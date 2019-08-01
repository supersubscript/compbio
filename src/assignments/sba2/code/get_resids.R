# library("Rpdb")
# pdb.files = list.files("../data", pattern = ".pdb")
# file = pdb.files[3]
# data = read.pdb(file)
# 
# ligand.name.1 = "IXO"
# ligand.name.2 = "CU2"
# 
# # is.ixo = which(data$atoms$resname == "IXO")
# amino.acids = c(
#   "GLY",
#   "ALA",
#   "VAL",
#   "LEU",
#   "ILE",
#   "PRO",
#   "PHE",
#   "TYR",
#   "TRP",
#   "SER",
#   "THR",
#   "CYS",
#   "MET",
#   "ASN",
#   "GLN",
#   "LYS",
#   "ARG",
#   "GLU"
# )
# 
# is.ixo = data$atoms$resname == "IXO"
# is.protein = data$atoms$resname %in% amino.acids
# dist = distances(data, is.ixo, is.protein)
# dist = norm(dist)
# 
# which()
# 

#
# data$atoms$resid
# # 
# in.file = "../data/3uon_within_5.dat"
# data = read.table(in.file)
# data = as.character(unique(data[,1]))
# write(data, file = paste0("", in.file, ".unique"), sep = "\n")
# # 
# # 3u
# CYS429 TRP400 TYR426 TYR430 ASP103 TYR104 TRP155 PHE181 THR190 THR187 ALA191 VAL407 ASN404 PHE195 ASN108 SER107
# # 4mqs
# PHE195 ALA194 VAl111 ASN108 TRP155 TYR104 ASP103 ASN404 TYR403 TRP400 CYS429 TYR426 TYR430  ASP103 ASN108 SER107 TRP400
# 
# # 4mqt_within_5_IXO.dat
# PHE195 ALA194 VAL111 ASN404 TRP400 TYR403 SER107 CYS429 TYR426 TYR430 ASP103 TYR104 TRP155 ASN108
# 
