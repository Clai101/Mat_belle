import os
os.system("mkdir root_mc")
os.system("mkdir hist_mc")
os.system("mkdir results_mc")
for root, dirs, files in os.walk("."):  
    for filename in files:
        if ".hist" in filename:
            os.system(f"h2root {filename}")
os.system(r"mv *.hist hist_mc")
os.system(r"mv *.root root_mc")

