library(abind) #load for binding 3D arrays when adding new segregating sites.

gen.multicopy.gcvar = function (pop, ne.n, cpndist.n = NULL, gc.rate = 0, gc.meanleng = 100) {
 #one generation of evolution: reproduction (with selection), reciprocal recombination, IGC, mutation, and monomorphic-site cleanup; returns an updated pop list.
 ne.o = (dim(pop$sites)[1])/2 #infer current Ne from the number of haplotypes (rows) in sites.
 ptnr = c(ne.o+(1:ne.o),1:ne.o) #fixed 1:1 haplotype pairing (1 to N+1, 2 to N+2, …) used for recombination template selection.
 loci = dim(pop$sites)[2]
 
 if (is.null(cpndist.n)) { #If new copy number distribution is not specified, reproduce with default model
  par1 = sample(1:(2*ne.o),2*ne.n,replace=TRUE,prob=c(pop$fit,pop$fit)) #parent sampling to form 2N gametes; selection enters via pop$fit.
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
 
 newsites = pop$sites[par1,,] #build offspring haplotypes by copying parental blocks (before recombination). 
 newcpn = pop$cpn[par1] #offspring copy numbers inherit from parents (may change after recombination).

 #Recombinations
 #A problem is the recom changes copy number, and it is impossible to set copy number dist AFTER recombination. This can only be remedied by checking whether copies are lost afterwards.
 #Note that pop$rec[1] is always assumed to be 0.
 
 rec.help = runif(2*ne.n,0,1) #per-haplotype uniform draws reused for all copy indices this generation.
 for (j in 2:dim(newsites)[3]) { #loop over copies 2 to k - reciprocal recombination: for haplotypes with rec.help < pop$rec[j], replace copy j with the paired partner’s copy j.
  rec.temp = which(rec.help < pop$rec[j])
  newsites[rec.temp,,j] = pop$sites[ptnr[par1[rec.temp]],,j]
 }
 for (i in which(rec.help < pop$rec[dim(newsites)[3]])) { #update newcpn - when recombination affects the last copy index, recompute copy number from the presence flags in column 1.
  newcpn[i] = sum(newsites[i,1,])
 }
 
 #Gene conversions
 gc.rate.i = gc.rate*(newcpn*(newcpn-1)) #scale IGC rate by ordered copy pairs per haplotype.
 gc.help = runif(2*ne.n,0,1)
 for (i in which(gc.help<gc.rate.i)) {
  gc.leng = rgeom(1,1/gc.meanleng) #draw IGC track length (geometrically distributed).
  if (gc.leng >= sleng) {gc.leng = sleng-1}
  gc.site = sample(sleng-gc.leng,1)
  temp = ((pop$pos >= gc.site)&(pop$pos < gc.site+gc.leng)) #sites that are within the GC track
  if (sum(temp) > 0) { #choose donor/recipient and overwrite track - interlocus gene conversion within a haplotype.
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
 if (sum(temp) > 0) { #append new segregating sites at the end of the loci array and assign unique positions in [1, sleng].
  newsites = abind(newsites,array(0,dim=c(2*ne.n,sum(temp),dim(newsites)[3])),along=2)
  for (i in 1:(2*ne.n)) {
   for (j in 1:dim(newsites)[3]) {
   if (temp[i,j]) {
    ct = ct + 1
	newsites[i,(loci+ct),j] = 1
   }
  }}
  
  for (i in 1:ct) { #draw unique genomic positions (1..sleng) for each new locus; avoid duplicates.
   mpos = sample(sleng,1)
   while (sum(newpos==mpos)>0) {mpos = sample(sleng,1)}
   newpos = c(newpos,mpos)
  }
 }

 #cleanup: compute per-locus allele frequency across copies, preserve two sentinel columns, and drop monomorphic sites.
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
 #function return - updated population with recalculated positions and copy numbers; fitness set to neutral (selection reapplied externally each gen).
 return(list(sites=newsites, pos=newpos, fit=newfit, cpn=newcpn, rec=pop$rec))

}

gen.write.multicopy.2 = function (pop, fname) {
 #write full population: one file per copy (haplotypes as g<0/1...>), SNP positions, and haplotype copy numbers.
 for (j in 1:dim(pop$sites)[3]) {
  out.sites = apply(cbind("g",pop$sites[,,j]),1,paste,collapse="")
  write(out.sites,paste(fname,"_copy_",j,".txt",sep=""),ncolumns=1)
 }
 write(pop$pos,paste(fname,"_pos.txt",sep=""),ncolumns=1)
 write(pop$cpn,paste(fname,"_copynumber.txt",sep=""),ncolumns=1)
}

gen.analysis.2copy.neutral = function (pop) {
 #This function analyzes the population and outputs key statistics, including frequencie of new copies, number of sites variable in each copy and in total, theta-pi and theta-w, fst and divergence between copies.
 hlc2 = (pop$sites[,1,2] == 1) #hlc = haplotype list containing copy 2

 freqs = matrix(0,length(pop$pos),2) #freqs: per-copy derived-allele counts at each locus; freqs.total collapses across all copies (used for “overall” summaries where the sample size is “all present copies”). 
 freqs[,1] = apply(pop$sites[,,1],2,sum)
 freqs[,2] = apply(pop$sites[,,2],2,sum)
 freqs.total = apply(pop$sites,2,sum)
  
 varsites = matrix(0,length(pop$pos),2) #varsites: marks loci variable within each copy using the relevant copy-specific sample size (2Ne for copy 1, number of carriers for copy 2).
 varsites[,1] = (freqs[,1] > 0) & (freqs[,1] < (2*ne))
 varsites[,2] = (freqs[,2] > 0) & (freqs[,2] < sum(hlc2))
  
 ovec = rep(0,16)
 
   ovec[1:3] = c(rr,gen,sum(hlc2)) #outputs replicate id, generation, counts of carriers of copy 2.
  
   ovec[4] = length(pop$pos)-2 #site counts. 4: total loci minus 2 sentinels; 5:6: variable sites in copies 1–2; 7: loci variable in both copies. 
   ovec[5:6] = apply(varsites[c(-1,-2),],2,sum)
   ovec[7] = sum(varsites[c(-1,-2),1]*varsites[c(-1,-2),2])
  
   ovec[8] = sum(2*freqs.total[c(-1,-2)]*(freqs.total[1]-freqs.total[c(-1,-2)])/(freqs.total[1]^2)) #8:9: overall theta_pi and theta_w using the total sample size (number of present copies, taken from the sentinel column via freqs.total[1]);
   ovec[9] = (length(pop$pos)-2)/sum(1/(1:freqs.total[1]))
   ovec[10] = sum(2*freqs[c(-1,-2),1]*((2*ne)-freqs[c(-1,-2),1])/((2*ne)^2)) #10:13: theta_pi and theta_w per copy (copy-specific sample sizes: 2Ne, sum(hlc2)).
   ovec[11] = sum(varsites[c(-1,-2),1])/sum(1/1:(2*ne))
   ovec[12] = sum(2*freqs[c(-1,-2),2]*(freqs[1,2]-freqs[c(-1,-2),2])/(freqs[1,2]^2))
   ovec[13] = sum(varsites[c(-1,-2),2])/sum(1/(1:freqs[1,2]))
  
   if (sum(hlc2) == 1) {ovec[14] = sum(pop$sites[hlc2,c(-1,-2),1] != pop$sites[hlc2,c(-1,-2),2])} #14: within-haplotype pairwise divergence (mean Hamming distance):
   else {ovec[14] = mean(apply((pop$sites[hlc2,c(-1,-2),1] != pop$sites[hlc2,c(-1,-2),2]),1,sum))}

  #FST (copy 1 vs copy 2): compute HS and HT with appropriate sample sizes, drop loci with HT=0, then output sum of per-locus FST (ovec[15]) and mean per-locus FST (ovec[16]).
   hs = ((2*freqs[c(-1,-2),1]*((2*ne)-freqs[c(-1,-2),1])/(2*ne)) + (2*freqs[c(-1,-2),2]*(freqs[1,2]-freqs[c(-1,-2),2])/freqs[1,2]))/((2*ne)+freqs[1,2])
   ht = 2*(freqs[c(-1,-2),1]+freqs[c(-1,-2),2])*((2*ne)+freqs[1,2]-freqs[c(-1,-2),1]-freqs[c(-1,-2),2])/(((2*ne)+freqs[1,2])^2)
   fsttemp = (ht >0)
   hs = hs[fsttemp]
   ht = ht[fsttemp]
   ovec[15] = sum((ht-hs)/ht)
   ovec[16] = mean((ht-hs)/ht)

 return(ovec)
 
}




arg = commandArgs(TRUE) #parse command-line parameters.

sel = as.numeric(arg[1]) #Positive selection coefficient for each copy beyond the original
gcr = as.numeric(arg[2]) #Here we use per-nucleotide IGC rate
gcml = as.numeric(arg[3]) #Mean length of IGC track
rec.rate = as.numeric(arg[4]) #Reciprocal recombination rate
rr = arg[5] #Replicate number

outhead = paste("Round1.5/Round1.5_adapt_",sel,"_rec_",rec.rate,"_GC_",gcr,"_length_",gcml,"_rep_",rr,sep="") #canonical output prefix including all parameters and replicate for reproducible filenames.
#core constants — Ne (diploid), sequence length, mu per bp, and mu per copy (mu.t = mu.bp * sleng)
ne = 2000
sleng = 20000
mu.bp = 5e-7
mu.t = mu.bp*sleng

#Read pop: read initial haplotypes - expect 2*ne lines of g + binary string; load into copy 1 (copy 2 empty).

inpop = read.table(paste("Neutral_Initial/Neutral_inieq_rep_",rr,".txt",sep=""),header=FALSE)
nn = dim(inpop)[1]
ini.loci = nchar(as.character(inpop$V1[1])) - 1
newsites.3d = array(0,dim=c(2*ne,(ini.loci+2),2))
 
for (i in 1:(2*ne)) {
 gntp.s = as.numeric(strsplit(as.character(inpop$V1[i]),"")[[1]][2:(ini.loci+1)])
 newsites.3d[i,,1] = c(1,1,gntp.s)
}
 
in.pos = c(-2,-1,sample(sleng,ini.loci)) #initialize positions - assign random positions to the initial ini.loci sites; sentinel positions are -2, -1.
 
pop.x = list(sites=newsites.3d,pos=in.pos,fit=rep(1,ne),cpn=rep(1,2*ne),rec=c(0,rec.rate)) #initial pop.x - set recombination vector as c(0, rec.rate) for the 2-copy stage (copy 1 never recombines).
#Read pop end



sw = 0

while (sw == 0) { #main loop with reset logic - evolve up to 200,000 generations; if copy 2 is lost (diploid copy sum == 2N) write a fail record and restart.
 sw = 1
 pop.run = pop.x
 
 #Here the first haplotype with second copy is randomized. This is to avoid having a low-fitness haplotype happen to get the copy which would make it hard to fix. Same for third.
 randhap = sample(2*ne,1)
 pop.run$sites[randhap,,2] = pop.run$sites[randhap,,1]
 pop.run$cpn[randhap] = 2
 copy.freq = NULL
 #write header - column names for the *_stats.txt summary file.
 write("Rep Gen Freq_2 Sites Sites_1 Sites_2 Sites_12 Tpi_all Tw_all Tpi_1 Tw_1 Tpi_2 Tw_2 Dist_12 sFst_12 Fst_12",paste(outhead,"_stats.txt",sep=""),ncolumns=1)
 
 for (gen in 1:200000) {
  pop.run$fit = (1+sel)^(pop.run$cpn[1:ne]+pop.run$cpn[ne+(1:ne)]-2) #selection each generation - diploid fitness (1+sel)^{(copies_on_diploid − 2)}; neutral when sel=0.

  pop.run = gen.multicopy.gcvar(pop.run,ne,gc.rate=gcr/gcml,gc.meanleng=gcml) #evolution step - call one generation with IGC rate specified per nucleotide: gc.rate = gcr/gcml converts to a per-track initiation rate.
  
  if (sum(pop.run$cpn) == (2*ne)) { #Return to initial state if the second copy is lost
   sw = 0
   write(gen,paste(outhead,"_failrecord_copy2.txt",sep=""),append=TRUE)
   break
  }
  
  if (gen %% 100 == 0) { #Analyze the population every 100 generations and output the statistics
   outvec = gen.analysis.2copy.neutral(pop.run)
   write(outvec,paste(outhead,"_stats.txt",sep=""),ncolumns=length(outvec),append=TRUE)
  }

 }
 if (sw == 1) {gen.write.multicopy.2(pop.run,paste(outhead,"_gen_200000",sep=""))} #final snapshot - write full population at gen 200,000 if the duplicate survives.
 

}



