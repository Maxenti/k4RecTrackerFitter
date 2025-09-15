total=1; perjob=1; seed0=4242
with open("joblist.txt","w") as f:
  for j in range((total+perjob-1)//perjob):
    ev=min(perjob,total-j*perjob)
    skip=j*perjob
    seed=seed0+j
    f.write(f"{j} {ev} {skip} {seed}\n")
print("Wrote joblist.txt")