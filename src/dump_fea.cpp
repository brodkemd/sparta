/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "dump_fea.h"

using namespace SPARTA_NS;

// customize by adding keyword

enum{INT,DOUBLE,BIGINT,STRING};        // same as Dump

enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain

/**
 * Execute system command and get STDOUT result.
 * Regular system() only gives back exit status, this gives back output as well.
 * @param command system command to execute
 * @return commandResult containing STDOUT (not stderr) output & exitstatus
 * of command. Empty if command failed (or has no output). If you want stderr,
 * use shell redirection (2&>1).
 */
static CommandResult EXEC(const std::string &command) {
    int exitcode = 0;
    std::array<char, 1048576> buffer {};
    std::string result;

    FILE *pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        throw std::runtime_error("popen() failed!");
    }
    try {
        std::size_t bytesread;
        while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
            result += std::string(buffer.data(), bytesread);
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    exitcode = WEXITSTATUS(pclose(pipe));
    return CommandResult{result, exitcode};
}

/**
 * Constructor
 */
DumpFea::DumpFea(SPARTA *sparta, int narg, char **arg) : DumpSurf(sparta, narg, arg) {
    if (narg == 5) error->all(FLERR,"No dump surf attributes specified");

    clearstep = 1;
    buffer_allow = 1;
    buffer_flag = 1;

    dimension = domain->dimension;

    int igroup = surf->find_group(arg[2]);
    if (igroup < 0) error->all(FLERR,"Dump surf group ID does not exist");
    groupbit = surf->bitmask[igroup];

    nevery = atoi(arg[3]);

    // expand args if any have wildcard character "*"
    // ok to include trailing optional args,
    //   so long as they do not have "*" between square brackets
    // nfield may be shrunk below if extra optional args exist

    int expand = 0;
    char **earg;
    int nargnew = nfield = input->expand_args(narg-5,&arg[5],1,earg);

    if (earg != &arg[5]) expand = 1;

    // allocate field vectors

    pack_choice = new FnPtrPack[nfield];
    vtype = new int[nfield];
    field2index = new int[nfield];
    argindex = new int[nfield];

    // custom props, computes, fixes, variables which the dump accesses

    ncustom = 0;
    id_custom = NULL;
    custom = NULL;

    ncompute = 0;
    id_compute = NULL;
    compute = NULL;

    nfix = 0;
    id_fix = NULL;
    fix = NULL;

    nvariable = 0;
    id_variable = NULL;
    variable = NULL;
    vbuf = NULL;

    // process attributes
    // ioptional = start of additional optional args in expanded args

    int ioptional = parse_fields(nfield,earg);

    if (ioptional < nfield)
        error->all(FLERR,"Invalid attribute in dump surf command");

    // noptional = # of optional args
    // reset nfield to subtract off optional args
    // reset ioptional to what it would be in original arg list

    int noptional = nfield - ioptional;
    nfield -= noptional;
    size_one = nfield;
    ioptional = narg - noptional;

    // setup format strings

    vformat = new char*[nfield];

    format_default = new char[4*nfield+1];
    format_default[0] = '\0';

    for (int i = 0; i < nfield; i++) {
        if (vtype[i] == INT) strcat(format_default,"%d ");
        else if (vtype[i] == DOUBLE) strcat(format_default,"%g ");
        else if (vtype[i] == BIGINT) strcat(format_default,BIGINT_FORMAT " ");
        vformat[i] = NULL;
    }

    format_column_user = new char*[size_one];
    for (int i = 0; i < size_one; i++) format_column_user[i] = NULL;

    // setup column string

    int n = 0;
    for (int iarg = 0; iarg < nfield; iarg++) n += strlen(earg[iarg]) + 2;
    columns = new char[n];
    columns[0] = '\0';
    for (int iarg = 0; iarg < nfield; iarg++) {
        strcat(columns,earg[iarg]);
        strcat(columns," ");
    }

    // if wildcard expansion occurred, free earg memory from expand_args()
    // wait to do this until after column string is created

    if (expand) {
        for (int i = 0; i < nargnew; i++) delete [] earg[i];
        memory->sfree(earg);
    }

    // trigger setup of list of owned surf elements belonging to surf group

    firstflag = 1;
    cglobal = clocal = NULL;
    buflocal = NULL;
}

/* ---------------------------------------------------------------------- */

/* Destructor */
DumpFea::~DumpFea() { 
    // syncing the file to update the handle
    fsync(fileno(this->fp));

    // clearing the string command to dump memory
    this->command.clear();
}

/* ---------------------------------------------------------------------- */

/**
 * overrides base class modify command 
 */ 
void DumpFea::modify_params(int narg, char** arg) {
    for (int i = 0; i < narg; i++) {
        // gets the command arg from the input
        if (arg[i] == (char*)"command")
            this->command = std::string(arg[i+1]);
    }
    return;
}

/**
 * writes to the dump file
 */
void DumpFea::write() {
    // opening the file for writing
    fp = fopen(this->filename, "w");

    // making sure it is open
    if (fp == NULL) error->all(FLERR,"Cannot open dump file");

    // simulation box bounds
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];

    // nme = # of dump lines this proc will contribute to dump
    nme = count();

    // ntotal = total # of dump lines in snapshot
    // nmax = max # of dump lines on any proc
    bigint bnme = nme;
    MPI_Allreduce(&bnme,&ntotal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

    int nmax;
    if (multiproc != nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
    else nmax = nme;

    // write timestep header
    // for multiproc,
    //   nheader = # of lines in this file via Allreduce on clustercomm
    bigint nheader = ntotal;
    if (multiproc)
        MPI_Allreduce(&bnme,&nheader,1,MPI_SPARTA_BIGINT,MPI_SUM,clustercomm);

    if (filewriter) write_header(nheader);

    // insure buf is sized for packing and communicating
    // use nmax to insure filewriter proc can receive info from others
    // limit nmax*size_one to int since used as arg in MPI calls
    if (nmax > maxbuf) {
        if ((bigint) nmax * size_one > MAXSMALLINT)
            error->all(FLERR,"Too much per-proc info for dump");
        maxbuf = nmax;
        memory->destroy(buf);
        memory->create(buf,maxbuf*size_one,"dump:buf");
    }

    // pack my data into buf
    pack();

    // if buffering, convert doubles into strings
    // insure sbuf is sized for communicating
    // cannot buffer if output is to binary file
    if (buffer_flag && !binary) {
        nsme = convert_string(nme,buf);
        int nsmin,nsmax;
        MPI_Allreduce(&nsme,&nsmin,1,MPI_INT,MPI_MIN,world);
        if (nsmin < 0) error->all(FLERR,"Too much buffered per-proc info for dump");
        if (multiproc != nprocs)
            MPI_Allreduce(&nsme,&nsmax,1,MPI_INT,MPI_MAX,world);
        else nsmax = nsme;
        if (nsmax > maxsbuf) {
            maxsbuf = nsmax;
            memory->grow(sbuf,maxsbuf,"dump:sbuf");
        }
    }

    // filewriter = 1 = this proc writes to file
    // ping each proc in my cluster, receive its data, write data to file
    // else wait for ping from fileproc, send my data to fileproc
    int tmp,nlines,nchars;
    MPI_Status status;
    MPI_Request request;

    // comm and output buf of doubles
    if (buffer_flag == 0 || binary) {
        if (filewriter) {
            for (int iproc = 0; iproc < nclusterprocs; iproc++) {
                if (iproc) {
                    MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,me+iproc,0,world,&request);
                    MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
                    MPI_Wait(&request,&status);
                    MPI_Get_count(&status,MPI_DOUBLE,&nlines);
                    nlines /= size_one;
                } else nlines = nme;

                write_data(nlines,buf);
            }
            // if (flush_flag) fflush(fp);

        } else {
            MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
            MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
        }

    // comm and output sbuf = one big string of formatted values per proc
    } else {
        if (filewriter) {
            for (int iproc = 0; iproc < nclusterprocs; iproc++) {
                if (iproc) {
                    MPI_Irecv(sbuf,maxsbuf,MPI_CHAR,me+iproc,0,world,&request);
                    MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
                    MPI_Wait(&request,&status);
                    MPI_Get_count(&status,MPI_CHAR,&nchars);
                } else nchars = nsme;

                write_data(nchars,(double *) sbuf);
            }
            // if (flush_flag) fflush(fp);

        } else {
            MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
            MPI_Rsend(sbuf,nsme,MPI_CHAR,fileproc,0,world);
        }
    }

    // closing the file
    fclose(this->fp);

    // if a command was provided
    if (this->command.length() > 0) {
        // running the command
        CommandResult command_result = EXEC(this->command);

        // if the command did not succeed
        if (command_result.exitstatus) {
            // writing to the logfile and erroring out
            fprintf(logfile, command_result.output.c_str());
            error->all(FLERR, "fix fea failed, see sparta log file");
        }
        std::cout << "done running command\n";
    }

}