library(abind)

gen.multicopy.gcvar = function (pop, ne.n, cpndist.n = NULL, gc.rate = 0, gc.meanleng = 100) {
 #This function evolves the population for one generation, doing reproduction based on fitness, then reciprocal recombination, IGC, mutation, and finally clean up the sites that are no longer polymorphic.
 ne.o = (dim(pop$sites)[1])/2
 ptnr = c(ne.o+(1:ne.o),1:ne.o)
 loci = dim(pop$sites)[2]
 
 if (is.null(cpndist.n)) { #If new copy number distribution is not specified, reproduce with default model
  par1 = sample(1:(2*ne.o),2*ne.n,replace=TRUE,prob=c(pop$fit,pop$fit)) #Parentage sampling: from haplotype to gamete.
 }
 else { #for example, cpndist.n = c(2000,1000,1000) is 2000 haplotypes with 1 copy, 1000 with 2 copies, 1000 with 3 copies
  par1 = NULL
  if (sum(cpndist.n) != 2*ne.n) {break}
  for (i in 1:length(cpndist.n)) {
   if ((cpndist.n[i]>0)&&(sum(pop$cpn==i)>0)) {
    if (sum(pop$cpn==i) == 1) {par1 = c(par1,rep(which(pop$cpn==i),cpndist.n[i]))}
    else {par1 = c(par1,sample(which(pop$cpn==i),cpndist.n[i],replace=TRUE,prob=c(pop$fit,pop$fit)[which(pop$cpn==i)]))}
   }
  }
  if (length(par1) != 2*ne.n) {break}
  par1 = sample(par1)
 }
 
 newsites = pop$sites[par1,,]
 newcpn = pop$cpn[par1]

 #Recombinations
 #A problem is the recom changes copy number, and it is impossible to set copy number dist AFTER recombination. This can only be remedied by checking whether copies are lost afterwards.
 #Note that pop$rec[1] is always assumed to be 0.
 
 rec.help = runif(2*ne.n,0,1)
 for (j in 2:dim(newsites)[3]) {
  rec.temp = which(rec.help < pop$rec[j])
  newsites[rec.temp,,j] = pop$sites[ptnr[par1[rec.temp]],,j]
 }
 for (i in which(rec.help < pop$rec[dim(newsites)[3]])) {
  newcpn[i] = sum(newsites[i,1,])
 }
 
 #Gene conversions
 gc.rate.i = gc.rate*(newcpn*(newcpn-1)) #GC rate corrected for copy number and direction.
 gc.help = runif(2*ne.n,0,1)
 for (i in which(gc.help<gc.rate.i)) {
  gc.leng = rgeom(1,1/gc.meanleng)
  if (gc.leng >= sleng) {gc.leng = sleng-1}
  gc.site = sample(sleng-gc.leng,1)
  temp = ((pop$pos >= gc.site)&(pop$pos < gc.site+gc.leng)) #sites that are within the GC track
  if (sum(temp) > 0) {
   gc.cp = sample(which(newsites[i,1,]==1),2) #randomly choose GC donor and recepient
   newsites[i,temp,gc.cp[2]] = newsites[i,temp,gc.cp[1]]
   
  }
 }
 
 #Mutations
 temp = (runif(2*ne.n*dim(newsites)[3],0,1) < mu.t) #Which COPIES have mutation? Non-existent copies are included, will be removed later
 temp = matrix(temp,2*ne.n,dim(newsites)[3]) #In matrix form. Rows -> first dimension in sites, columns -> third dimension
 temp = temp*newsites[,1,] #All mutations in non-existent copies are multiplied with zero
 ct = 0 #Counting number of individuals already mutated, to know which locus to mutate.
 newpos = pop$pos
 if (sum(temp) > 0) {
  newsites = abind(newsites,array(0,dim=c(2*ne.n,sum(temp),dim(newsites)[3])),along=2)
  for (i in 1:(2*ne.n)) {
   for (j in 1:dim(newsites)[3]) {
   if (temp[i,j]) {
    ct = ct + 1
	newsites[i,(loci+ct),j] = 1
   }
  }}
  
  for (i in 1:ct) {
   mpos = sample(sleng,1)
   while (sum(newpos==mpos)>0) {mpos = sample(sleng,1)}
   newpos = c(newpos,mpos)
  }
 }

 #afq cleanup
 temp = NULL
 for (i in 1:dim(newsites)[3]) {temp = rbind(temp,newsites[,,i])}
 temp = temp[(temp[,1]+temp[,2]>0),] #Remove lines corresponding to non-existent copies
 afq = apply(temp,2,mean)
 afq[1] = 0.5
 afq[2] = 0.5
 if (sum(afq) != 0 & sum(afq>0) != 1) {
  newsites = newsites[,(afq*(afq-1) != 0),]
  newpos = newpos[(afq*(afq-1) != 0)]
 }
 loci = dim(newsites)[2]


 #Calculate fitness (none at this time)
 newfit = rep(1,ne.n)
 
 return(list(sites=newsites, pos=newpos, fit=newfit, cpn=newcpn, rec=pop$rec))

}

gen.write.multicopy.2 = function (pop, fname) {
 #This function writes out the entire population, including alleles on each snp for each haplotype, positions of snps, and copy number for each haplotype.
 for (j in 1:dim(pop$sites)[3]) {
  out.sites = apply(cbind("g",pop$sites[,,j]),1,paste,collapse="")
  write(out.sites,paste(fname,"_copy_",j,".txt",sep=""),ncolumns=1)
 }
 write(pop$pos,paste(fname,"_pos.txt",sep=""),ncolumns=1)
 write(pop$cpn,paste(fname,"_copynumber.txt",sep=""),ncolumns=1)
}

gen.analysis.3copy.neutral = function (pop) {
 #This function analyzes the population and outputs key statistics, including frequencie of new copies, number of sites variable in each copy and in total, theta-pi and theta-w, fst and divergence among the three copies.

 if (max(pop$cpn) == 3) {sw3 = 1}
 else {sw3 = 0}

 hlc2 = (pop$sites[,1,2] == 1) #hlc = haplotype list containing (copy 2 or 3)
 hlc3 = rep(FALSE,length(hlc2))
 hlc23 = rep(FALSE,length(hlc2))
 if (sw3) {
  hlc3 = (pop$sites[,1,3] == 1)
  hlc23 = (hlc2&hlc3)
 }
 
 freqs = matrix(0,length(pop$pos),3) #Frequency by count, each copy separately
 freqs[,1] = apply(pop$sites[,,1],2,sum)
 freqs[,2] = apply(pop$sites[,,2],2,sum)
 if (sw3) {freqs[,3] = apply(pop$sites[,,3],2,sum)}
 freqs.total = apply(pop$sites,2,sum)

 varsites = matrix(0,length(pop$pos),3) #Boolean for whether a site is variable within copy
 varsites[,1] = (freqs[,1] > 0) & (freqs[,1] < 4000)
 varsites[,2] = (freqs[,2] > 0) & (freqs[,2] < sum(hlc2))
 if (sw3) {varsites[,3] = (freqs[,3] > 0) & (freqs[,3] < sum(hlc3))}
 
 ovec = rep(0,32)

  ovec[1:5] = c(rr,gen,sum(hlc2),sum(hlc3),sum(hlc23))
  ovec[6] = length(pop$pos)-2
  ovec[7:9] = apply(varsites,2,sum)
  ovec[10] = sum(varsites[,1]*varsites[,2])
  ovec[11] = sum(varsites[,1]*varsites[,3])
  ovec[12] = sum(varsites[,2]*varsites[,3])
  ovec[13] = sum(varsites[,1]*varsites[,2]*varsites[,3]) #col6~13: number of variable sites, total, each copy, copy shared
  
  ovec[14] = sum(2*freqs.total*(freqs.total[1]-freqs.total)/(freqs.total[1]^2))
  ovec[15] = (length(pop$pos)-2)/sum(1/(1:freqs.total[1]))
  ovec[16] = sum(2*freqs[,1]*(4000-freqs[,1])/(4000^2))
  ovec[17] = sum(varsites[,1])/8.87139
  ovec[18] = sum(2*freqs[,2]*(freqs[1,2]-freqs[,2])/(freqs[1,2]^2))
  ovec[19] = sum(varsites[,2])/sum(1/(1:freqs[1,2]))
  if (sw3) {
   ovec[20] = sum(2*freqs[,3]*(freqs[1,3]-freqs[,3])/(freqs[1,3]^2))
   ovec[21] = sum(varsites[,3])/sum(1/(1:freqs[1,3])) #col14~21: theta_pi and theta_w for the whole dataset and each copy
  }

  if (sum(hlc2) == 1) {ovec[22] = sum(pop$sites[hlc2,,1] != pop$sites[hlc2,,2])}
  else {ovec[22] = mean(apply((pop$sites[hlc2,,1] != pop$sites[hlc2,,2]),1,sum))}
   
  if (sw3) {
   if (sum(hlc3) == 1) {ovec[23] = sum(pop$sites[hlc3,,1] != pop$sites[hlc3,,3])}
   else {ovec[23] = mean(apply((pop$sites[hlc3,,1] != pop$sites[hlc3,,3]),1,sum))}
   if (sum(hlc23) == 1) {ovec[24] = sum(pop$sites[hlc23,,2] != pop$sites[hlc23,,3])}
   else {ovec[24] = mean(apply((pop$sites[hlc23,,2] != pop$sites[hlc23,,3]),1,sum))}
  }
  
  hs = ((2*freqs[,1]*(4000-freqs[,1])/4000) + (2*freqs[,2]*(freqs[1,2]-freqs[,2])/freqs[1,2]))/(4000+freqs[1,2])
  ht = 2*(freqs[,1]+freqs[,2])*(4000+freqs[1,2]-freqs[,1]-freqs[,2])/((4000+freqs[1,2])^2)
  temp = (ht >0)
  hs = hs[temp]
  ht = ht[temp]
  ovec[25] = sum((ht-hs)/ht)
  ovec[26] = mean((ht-hs)/ht) #col25~26: Fst between copy 1 and 2
  
  if (sw3) {
   hs = ((2*freqs[,1]*(4000-freqs[,1])/4000) + (2*freqs[,3]*(freqs[1,3]-freqs[,3])/freqs[1,3]))/(4000+freqs[1,3])
   ht = 2*(freqs[,1]+freqs[,3])*(4000+freqs[1,3]-freqs[,1]-freqs[,3])/((4000+freqs[1,3])^2)
   temp = (ht >0)
   hs = hs[temp]
   ht = ht[temp]
   ovec[27] = sum((ht-hs)/ht)
   ovec[28] = mean((ht-hs)/ht) #col27~28: Fst between copy 1 and 3
  
   hs = ((2*freqs[,2]*(freqs[1,2]-freqs[,2])/freqs[1,2]) + (2*freqs[,3]*(freqs[1,3]-freqs[,3])/freqs[1,3]))/(freqs[1,2]+freqs[1,3])
   ht = 2*(freqs[,2]+freqs[,3])*(freqs[1,2]+freqs[1,3]-freqs[,2]-freqs[,3])/((freqs[1,2]+freqs[1,3])^2)
   temp = (ht >0)
   hs = hs[temp]
   ht = ht[temp]
   ovec[29] = sum((ht-hs)/ht)
   ovec[30] = mean((ht-hs)/ht) #col29~30: Fst between copy 2 and 3
   
   hs = ((2*freqs[,1]*(4000-freqs[,1])/4000) + (2*freqs[,2]*(freqs[1,2]-freqs[,2])/freqs[1,2]) + (2*freqs[,3]*(freqs[1,3]-freqs[,3])/freqs[1,3]))/(4000+freqs[1,2]+freqs[1,3])
   ht = 2*freqs.total*(freqs.total[1]-freqs.total)/(freqs.total[1]^2)
   temp = (ht >0)
   hs = hs[temp]
   ht = ht[temp]
   ovec[31] = sum((ht-hs)/ht)
   ovec[32] = mean((ht-hs)/ht) #col31~32: Fst of all copies
  }
  else {
   ovec[31:32] = ovec[25:26]
  }
  
  return(ovec)

}

arg = commandArgs(TRUE)

sel = as.numeric(arg[1]) #Positive selection coefficient for each copy beyond the original
gcr = as.numeric(arg[2]) #Here we use per-nucleotide IGC rate
gcml = as.numeric(arg[3]) #Mean length of IGC track
dup.time = as.numeric(arg[4]) #Time when the third copy is created by a second duplication event
rec.rate = as.numeric(arg[5]) #Reciprocal recombination rate
rr = arg[6] #Replicate number

outhead = paste("Round2.5/Round2.5_adapt_",sel,"_rec_",rec.rate,"_GC_",gcr,"_length_",gcml,"_timeitv_",dup.time,"_rep_",rr,sep="")

ne = 2000
sleng = 20000
mu.bp = 5e-7
mu.t = mu.bp*sleng

#Read pop

inpop = read.table(paste("Neutral_Initial/Neutral_inieq_rep_",rr,".txt",sep=""),header=FALSE)
nn = dim(inpop)[1]
ini.loci = nchar(as.character(inpop$V1[1])) - 1
newsites.3d = array(0,dim=c(2*ne,(ini.loci+2),2))
 
for (i in 1:(2*ne)) {
 gntp.s = as.numeric(strsplit(as.character(inpop$V1[i]),"")[[1]][2:(ini.loci+1)])
 newsites.3d[i,,1] = c(1,1,gntp.s)
}
 
in.pos = c(-2,-1,sample(sleng,ini.loci))
 
pop.x = list(sites=newsites.3d,pos=in.pos,fit=rep(1,ne),cpn=rep(1,2*ne),rec=c(0,rec.rate))
#Read pop end



sw = 0

while (sw == 0) {
 sw = 1
 pop.run = pop.x
 
 #Here the first haplotype with second copy is randomized. This is to avoid having a low-fitness haplotype happen to get the copy which would make it hard to fix. Same for third.
 randhap = sample(2*ne,1)
 pop.run$sites[randhap,,2] = pop.run$sites[randhap,,1]
 pop.run$cpn[randhap] = 2
 copy.freq = NULL

 write("Rep Gen Freq_2 Freq_3 Freq_23 Sites Sites_1 Sites_2 Sites_3 Sites_12 Sites_13 Sites_23 Sites_123 Tpi_all Tw_all Tpi_1 Tw_1 Tpi_2 Tw_2 Tpi_3 Tw_3 Dist_12 Dist_13 Dist_23 sFst_12 Fst_12 sFst_13 Fst_13 sFst_23 Fst_23 sFst_all Fst_all",paste(outhead,"_stats.txt",sep=""),ncolumns=1)

 for (gen in 1:dup.time) {
  pop.run$fit = (1+sel)^(pop.run$cpn[1:ne]+pop.run$cpn[ne+(1:ne)]-2)
  pop.run = gen.multicopy.gcvar(pop.run,ne,gc.rate=gcr/gcml,gc.meanleng=gcml)
  
  if (sum(pop.run$cpn) == 2*ne) { #Return to initial state if the second copy is lost
   sw = 0
   write(gen,paste(outhead,"_failrecord_copy2.txt",sep=""),append=TRUE)
   break
  }
  
  if (gen %% 100 == 0) { #Analyze the population every 100 generations and output the statistics
   outvec = gen.analysis.3copy.neutral(pop.run)
   write(outvec,paste(outhead,"_stats.txt",sep=""),ncolumns=length(outvec),append=TRUE)
  }

  
 }

}


#Create the third copy

pop.x2 = pop.run

pop.x2$sites = abind(pop.x2$sites,matrix(0,2*ne,dim(pop.x2$sites)[2]),along=3)
pop.x2$rec = c(0,rec.rate,2*rec.rate)

sw = 0

while (sw == 0) {
 sw = 1
 pop.run = pop.x2

 randhap = sample(which(pop.run$sites[,1,2]==1),1)
 pop.run$sites[randhap,,3] = pop.run$sites[randhap,,2]
 pop.run$cpn[randhap] = 3

 
 for (gen in (dup.time+1):200000) {
  pop.run$fit = (1+sel)^(pop.run$cpn[1:ne]+pop.run$cpn[ne+(1:ne)]-2)
  pop.run = gen.multicopy.gcvar(pop.run,ne,gc.rate=gcr/gcml,gc.meanleng=gcml)
  
  if ((sum(pop.run$sites[,1,2])*sum(pop.run$sites[,1,3])) == 0) {#Return to the time of second duplication event, if either second or third copy is lost
   sw = 0
   write(gen,paste(outhead,"_failrecord_copy3.txt",sep=""),append=TRUE)
   break
  }
  
  if (gen %% 100 == 0) { #Analyze the population every 100 generations and output the statistics
   outvec = gen.analysis.3copy.neutral(pop.run)
   write(outvec,paste(outhead,"_stats.txt",sep=""),ncolumns=length(outvec),append=TRUE)
  }
  
 }

 if (sw == 1) {gen.write.multicopy.2(pop.run,paste(outhead,"_gen_200000",sep=""))}


}






