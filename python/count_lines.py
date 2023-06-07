import os

check_dir = os.path.sep + os.path.join(*__file__.split(os.path.sep)[:-2], "src")

count = 0
num_files = 0
for item in os.listdir(check_dir):
    if item.endswith(".h") or item.endswith(".cpp"):
        item = os.path.join(check_dir, item)
        num_files+=1
        print("counting:", item)
        with open(item, 'r') as f: count += len(f.readlines())

print(num_files, "files,", count, "lines")