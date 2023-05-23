from ._src import _checkType


class Species:
    def __init__(self,
            id:str,
            molweight:str|float, # molecular weight in amu (atomic mass units, e.g. 16 for oxygen)
            molmass:str|float, # molecular mass (mass units)
            rotational_dof:str|int, # rotational degrees of freedom (integer, unitless)
            inv_rotational_relaxation:str|float, # inverse rotational relaxtion number (unitless)
            vibrational_dof:str|int, # vibrational degrees of freedom (integer, unitless)
            inv_vibrational_relaxation:str|float, # inverse vibrational relaxation number
            vibrational_temp:str|float, # vibrational temperature (temperature units)
            species_wt:str|float, # species weight (unitless)
            charge:str|float, # multiple of electon charge (1 for a proton)
        ):
        
        for item in [id, molweight, molmass, rotational_dof, inv_rotational_relaxation, vibrational_dof, inv_vibrational_relaxation, vibrational_temp, species_wt, charge]:
            _checkType("input must be str, float, int", item, str, float, int)
        self.id:str = id
        self.molweight:str|float = molweight
        self.molmass:str|float = molmass
        self.rotational_dof:str|int = rotational_dof
        self.inv_rotational_relaxation:str|float = inv_rotational_relaxation
        self.vibrational_dof:str|int = vibrational_dof
        self.inv_vibrational_relaxation:str|float = inv_vibrational_relaxation
        self.vibrational_temp:str|float = vibrational_temp
        self.species_wt:str|float = species_wt
        self.charge:str|float = charge

    def __call__(self) -> tuple:
        return self.id, self.molweight, self.molmass, self.rotational_dof, self.inv_rotational_relaxation, self.vibrational_dof, self.inv_vibrational_relaxation, self.vibrational_temp, self.species_wt, self.charge
    
    def __str__(self) -> str:
        _s = ""
        for item in [self.id, self.molweight, self.molmass, self.rotational_dof, self.inv_rotational_relaxation, self.vibrational_dof, self.inv_vibrational_relaxation, self.vibrational_temp, self.species_wt, self.charge]:
            _s += str(item) + " "

        return _s.strip()

class VSS:
    def __init__(self,
        id:str,
        diameter:str|float, #  VHS or VSS diameter of particle (distance units)
        omega:str|float, # temperature-dependence of viscosity (unitless)
        tref:str|float, # reference temperature (temperature units)
        alpha:str|float # angular scattering parameter (unitless)
    ) -> None:
        
        for item in [diameter, omega, tref, alpha]:
            _checkType("input must be str, float, int", item, str, float, int)
        
        self.id:str = id
        self.diameter:str|float = diameter
        self.omega:str|float = omega
        self.tref:str|float = tref
        self.alpha:str|float = alpha

    def __call__(self) -> tuple:
        return self.id, self.diameter, self.omega, self.tref, self.alpha

    def __str__(self) -> str:
        _s = ""
        for item in [self.id, self.diameter, self.omega, self.tref, self.alpha]:
            _s += str(item) + " "

        return _s.strip()


class Molecule:
    def __init__(self, 
            id:str, # exp. for oxygen this is O2 
            molweight:str|float, # molecular weight in amu (atomic mass units, e.g. 16 for oxygen)
            molmass:str|float, # molecular mass (mass units)
            rotational_dof:str|int, # rotational degrees of freedom (integer, unitless)
            inv_rotational_relaxation:str|float, # inverse rotational relaxtion number (unitless)
            vibrational_dof:str|int, # vibrational degrees of freedom (integer, unitless)
            inv_vibrational_relaxation:str|float, # inverse vibrational relaxation number
            vibrational_temp:str|float, # vibrational temperature (temperature units)
            species_wt:str|float, # species weight (unitless)
            charge:str|float, # multiple of electon charge (1 for a proton)
            diameter:str|float, #  VHS or VSS diameter of particle (distance units)
            omega:str|float, # temperature-dependence of viscosity (unitless)
            tref:str|float, # reference temperature (temperature units)
            alpha:str|float # angular scattering parameter (unitless)
        ):

        if not (isinstance(id, str)):
            raise TypeError(f"id must be string")
        
        self.id:str = id
        self.species:Species = Species(
            id,
            molweight,
            molmass,
            rotational_dof,
            inv_rotational_relaxation,
            vibrational_dof,
            inv_vibrational_relaxation,
            vibrational_temp,
            species_wt,
            charge
        )

        self.vss:VSS = VSS(
            id,
            diameter,
            omega,
            tref,
            alpha
        )

if __name__ == "__main__":
    Species("s", 0, 0, 0, 0, 0, 0, 0, [])