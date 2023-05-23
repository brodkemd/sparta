import os, shutil

pkg_path = os.path.sep+os.path.join(*__file__.split(os.path.sep)[1:-1])
main_path = os.path.sep+os.path.join(*__file__.split(os.path.sep)[1:-2])

doc_src = os.path.join(main_path, "doc")
doc_dest = os.path.join(pkg_path, "src", "sparta", "doc")

try:
    shutil.rmtree(doc_dest)
except FileNotFoundError: pass

try:
    os.mkdir(doc_dest)
except FileExistsError: pass

for item in os.listdir(doc_src):
    entry_path = os.path.join(doc_src, item)
    
    if os.path.isdir(entry_path):
        dest_dir = entry_path.replace(doc_src, doc_dest)
        print("Making:", dest_dir)
        os.mkdir(dest_dir)

        for item_1 in os.listdir(entry_path):
            if item_1[item_1.rfind(".") + 1:] in ["tex", "txt"]:
                print("- Skipping:", item)
                continue
            entry_1_path = os.path.join(entry_path, item_1)
            dest_path = entry_1_path.replace(doc_src, doc_dest)
            # print("1:", entry_1_path)
            print("Copying:", item_1, " to ", dest_path)
            shutil.copy2(entry_1_path, dest_path)
            #exit()
    else:
        # print("0:", entry_path, end=" ")
        dest_path = entry_path.replace(doc_src, doc_dest)

        if item[item.rfind(".") + 1:] in ["tex", "txt"]:
            print("- Skipping:", item)
            continue
        
        print("Copying:", item)

        shutil.copy2(entry_path, dest_path)
    