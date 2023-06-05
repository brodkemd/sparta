import os

count = 0
num_files = 0
for item in os.listdir():
    if item.endswith(".h") or item.endswith(".cpp"):
        num_files+=1
        print("counting:", item)
        with open(item, 'r') as f: count += len(f.readlines())

print(num_files, "files,", count, "lines")