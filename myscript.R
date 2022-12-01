install.packages("ggplot2")
library(ggplot2)


x = 1.5
y = exp(x)
#print(y)


# Writing functions

myfunction = function(x,y)
{
z = x+y
return(z)
}

#Write a function for euclidian distance

Distance = function(x,y)
{
d = sqrt(sum(x-y)^2)
return(d)
}

#Fibonacci numbers
fib = function(n)
{
a = c(0,1)
for (i in 3:n)
{
a[i] = a[i-1] + a[i-2] 

}
return(a)
}

#Computer exercise 2 part 6 coverage
genome1.subset=genome1[1:1000,]
ref.subset=reference[1:1000]

sum_gen1 = apply(genome1.subset, 1, sum)
sum_gen1 = sum_gen1 - genome1.subset["Position"]
mean_gen1 = mean(sum_gen1[,1])

#For entire genome1
coverage1 = apply(genome1, 1, sum)
coverage1 = coverage1 - genome1["Position"]
mean_cov1 = mean(coverage1[,1])
max_cov1 = max(coverage1[,1])
min_cov1 = min(coverage1[,1])

int_cov1_1 = coverage1[,1][0:1000]
int_cov1_2 = coverage1[,1][1001:2000]



plot(int_cov1_1, main = "Coverage of first 1000 positions in the first genome",
     xlab = "Position", ylab = "Coverage",
     pch=20, frame = FALSE, col = "cornflowerblue", ylim = c(0,30))

plot(int_cov1_2, main = "Coverage of second 1000 positions in the first genome",
     xlab = "Position", ylab = "Coverage",
     pch=20, frame = FALSE, col = "cornflowerblue", ylim = c(0,30))

#Part 7: Error rate

genome1.length=nrow(genome1) # Calculate the length of genome1
# Allocate an empty vector for matches
matches=vector(length=genome1.length)
for (pos in 1:genome1.length)
{
  # Picks the number of reads that has the same nucleotide as
  # the reference
  matches[pos] = genome1[pos,reference[pos]]
}

mismatches = coverage1[,1] - matches

#The number of positions that at least have one read with a mismatch
d = 0
for (dx in 1:4641652)
{
  if (mismatches[dx] > 0)
  {d = d+1}
}

error = mismatches/coverage1
error2 = error[,1][1001:2000]

error3 = error2
error3[error3==0]<- NA


plot(error3,main = "Error rate over the second 1000 positions",
     xlab = "Position", ylab = "Error rate",
     pch=20, frame = FALSE, col = "blueviolet", ylim = c(0,0.1))


#Part8:A test for single nucleotide polymorphisms 
#perror is 0.01
binom_test = function(num_mis, num_cov)
{
 if (num_cov == 0)
 {return -1}
  else
  {b = binom.test(num_mis, num_cov, p = 0.01)}
  return(b)
}
#Small p value means mutation
#Part9:Screening the genome for SNPs

p_values = vector(length=1000)

for (ix in 1:1000)
{
  probabilities = binom_test(mismatches[ix],coverage1[,1][ix])
  p_values[ix] = probabilities[["p.value"]]
}
#Apply it to the whole genome one
p_tot_values = vector(length=genome1.length)

for (jx in 1:genome1.length)
{
  probability_tot = binom_test(mismatches[jx],coverage1[,1][jx])
  p_tot_values[jx] = probability_tot[["p.value"]]
  
}
#The positions of the lowest p values(mutations)
p_ordered = order(p_tot_values, decreasing = FALSE)


#Select the significant positions
k = 0
for (kx in 1:genome1.length)
{
  if (p_tot_values[kx] <= 0.05)
  {k = k+1}
}

#Checking the significant mutations

genome1[order(p_tot_values, decreasing = F)[1:10],]
reference[order(p_tot_values, decreasing = F)[1:10]]
sort(p_tot_values, decreasing = F)[1:100]

#Genome2
coverage2 = apply(genome2, 1, sum)
coverage2 = coverage2 - genome2["Position"]

genome2.length=nrow(genome2) # Calculate the length of genome2
# Allocate an empty vector for matches
matches2=vector(length=genome2.length)
for (pos in 1:genome2.length)
{
  # Picks the number of reads that has the same nucleotide as
  # the reference
  matches2[pos] = genome2[pos,reference[pos]]
}

mismatches2 = coverage2[,1] - matches2

p_tot_values2 = vector(length=genome2.length)
for (jx in 1:genome2.length)
{
  probability_tot2 = binom_test(mismatches2[jx],coverage2[,1][jx])
  p_tot_values2[jx] = probability_tot2[["p.value"]]
  
}

genome2[order(p_tot_values2, decreasing = F)[1:10],]
reference[order(p_tot_values2, decreasing = F)[1:10]]
sort(p_tot_values2, decreasing = F)[1:100]

#We found only one mutation in genome 2

#Genome3
coverage3 = apply(genome3, 1, sum)
coverage3 = coverage3 - genome3["Position"]

genome3.length=nrow(genome3) # Calculate the length of genome3
# Allocate an empty vector for matches
matches3=vector(length=genome3.length)
for (pos in 1:genome3.length)
{
  # Picks the number of reads that has the same nucleotide as
  # the reference
  matches3[pos] = genome3[pos,reference[pos]]
}

mismatches3 = coverage3[,1] - matches3

p_tot_values3 = vector(length=genome3.length)
for (jx in 1:genome3.length)
{
  probability_tot3 = binom_test(mismatches3[jx],coverage3[,1][jx])
  p_tot_values3[jx] = probability_tot3[["p.value"]]
  
}

genome3[order(p_tot_values3, decreasing = F)[1:10],]
reference[order(p_tot_values3, decreasing = F)[1:10]]
sort(p_tot_values3, decreasing = F)[1:100]
#We found three mutations
