import glob
files = glob.glob("*.gz")
import pandas

for i in xrange(len(files)):
   data = pandas.read_table(files[i],sep="\t",compression="gzip")
   if files[i].find("UKBB") == -1:
      data = data[data["EFFECT_ALLELE"] != "I"]
      data = data[data["EFFECT_ALLELE"] != "D"]
   data = data[data["IMPUTATION"] > 0.4]
   data = data[data["MAC"] > 5]
   a = data["cptid"].values
   chr = [a[j][0:a[j].index(":")] for j in xrange(len(a))]
   pos = [a[j][a[j].index(":")+1:] for j in xrange(len(a))]
   if files[i].find("UKBB") != -1:
      for j in xrange(len(pos)): 
         if pos[j].find(":") != -1:
            pos[j] = pos[j][0:pos[j].index(":")] 
   for j in xrange(len(chr)):
      if chr[j] == "X":
         chr[j] = "23"
   a1 = data["EFFECT_ALLELE"].values
   a2 = data["OTHER_ALLELE"].values
   uniq = []
   for j in xrange(len(a1)):
      if a1[j] == "A" and a2[j] == "C":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a1[j] + "_" + a2[j])
      elif a1[j] == "A" and a2[j] == "G":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a1[j] + "_" + a2[j])
      elif a1[j] == "A" and a2[j] == "T":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a1[j] + "_" + a2[j])
      elif a1[j] == "C" and a2[j] == "A":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a2[j] + "_" + a1[j])
      elif a1[j] == "C" and a2[j] == "G":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a1[j] + "_" + a2[j]) 
      elif a1[j] == "C" and a2[j] == "T":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a1[j] + "_" + a2[j])
      elif a1[j] == "G" and a2[j] == "A":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a2[j] + "_" + a1[j])
      elif a1[j] == "G" and a2[j] == "C":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a2[j] + "_" + a1[j])
      elif a1[j] == "G" and a2[j] == "T":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a1[j] + "_" + a2[j])
      elif a1[j] == "T" and a2[j] == "A":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a2[j] + "_" + a1[j])
      elif a1[j] == "T" and a2[j] == "C":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a2[j] + "_" + a1[j])
      elif a1[j] == "T" and a2[j] == "G":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a2[j] + "_" + a1[j])
      elif a1[j] == "I" and a2[j] == "D":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a2[j] + "_" + a1[j])
      elif a1[j] == "D" and a2[j] == "I":
         uniq.append(chr[j] + ":" + pos[j] + "_" + a1[j] + "_" + a2[j])
   data["Unique_ID"] = uniq
   chr = [int(float(chr[j])) for j in xrange(len(chr))]
   pos = [int(float(pos[j])) for j in xrange(len(pos))]
   data["CHR"] = chr
   data["BP"] = pos
   data["STRAND"] = "+"
   if "SNP" in data.columns:
      data = data.drop(axis=1,labels="SNP")
   data.to_csv(files[i][0:files[i].index(".gz")] + "_Preproc.gz",sep="\t",index=False,compression="gzip")
   print repr(i+1) + "/" + repr(len(files))