mutant_raw = read.csv("Moraes_blue.csv", header=FALSE)
wild_raw = read.csv("Moraes_red.csv", header=FALSE)

names(mutant_raw) = c("days", "mtDNA levels")
names(wild_raw) = c("days", "mtDNA levels")

par(mfrow=c(1,2))
plot(mutant_raw)
plot(wild_raw)
par(mfrow=c(1,1))

time = rep(c(0,6,15,22,30,45), 6)
type = c(rep("wild", 18), rep("mutant", 18))
bound = rep(c(rep("lower", 6), rep("pred", 6), rep("upper", 6)), 2)

Moraes_data = data.frame(time=time, exp.level=c(wild_raw[,2], mutant_raw[,2]), 
                         type=type, bound=bound)

write.csv(Moraes_data, file="Moraes_data.csv")

