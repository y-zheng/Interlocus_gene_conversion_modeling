gen.onegene.norec = function (pop, ne.n) { #single-copy generation: selection is applied outside; here we do parent sampling, mutations, and monomorphic-site cleanup.

 newsites = NULL
 ne.o = (dim(pop$sites)[1])/2
 ptnr = c(ne.o+(1:ne.o),1:ne.o)
 loci = dim(pop$sites)[2]
 par1 = sample(1:(2*ne.o),2*ne.n,replace=TRUE,prob=c(pop$fit,pop$fit)) #parent sampling to form 2N haploid gametes; selection acts via pop$fit.

 #recombination code is intentionally disabled for the single-gene stage.
 #rec.info = rbinom(2*ne.n,1,rec.rate)*runif(2*ne.n,0,1) #0 if no recombination; otherwise position of recombination
 #rec.sites = (t(matrix(pop$pos,loci,2*ne.n))>=matrix(rec.info,2*ne.n,loci))
 #newsites = pop$sites[par1,] * rec.sites + pop$sites[ptnr[par1],] * (1-rec.sites)
 newsites = pop$sites[par1,] #copy parental haplotypes (before mutation) to form the offspring matrix.
 
 #Mutations
 temp = (runif(2*ne.n,0,1) < mu.t) #Bernoulli per haplotype for a new mutation with rate mu.t = mu.bp * sleng.
 ct = 0 #running counter of new mutations this generation (to index new loci).
 newpos = pop$pos
 if (sum(temp) > 0) {
  newsites = cbind(newsites,matrix(0,2*ne.n,sum(temp))) #append new (all-zero) columns for each realized mutation.
  for (i in 1:(2*ne.n)) {
   if (temp[i]) {
    ct = ct + 1
	newsites[i,(loci+ct)] = 1 #set the focal haplotype’s new locus to the derived state (1).
   }
  }
  
  for (i in 1:ct) { #draw unique genomic positions (1..sleng) for each new locus; avoid duplicates.
   mpos = sample(sleng,1)
   while (sum(newpos==mpos)>0) {mpos = sample(sleng,1)}
   newpos = c(newpos,mpos)
  }
 }

 #start cleanup: mark fixed/lost loci for removal while preserving sentinels.
 afq = apply(newsites,2,mean) #per-locus allele frequency; set sentinels (cols 1–2) to 0.5 so they survive cleanup.
 afq[1] = 0.5
 afq[2] = 0.5
 if (sum(afq) != 0 & sum(afq>0) != 1) { #drop monomorphic loci (keeps only sites with 0 < afq < 1).
  newsites = newsites[,(afq*(afq-1) != 0)]
  newpos = newpos[(afq*(afq-1) != 0)]
 }
 loci = dim(newsites)[2]

 #intra-function fitness is neutral (1); selection is reapplied in the outer loop.
 newfit = rep(1,ne.n)
 
 return(list(sites=newsites, pos=newpos, fit=newfit)) #return updated population: new sites, positions, and neutral fitness.

}

gen.read.onegene = function (fname) { #read neutral haplotypes and add two sentinels.
 #input file has no sentinels; they’re added here to protect early columns from cleanup.
 inpop = read.table(fname,header=FALSE)
 nn = dim(inpop)[1]
 loci = nchar(as.character(inpop$V1[1])) - 1 #number of segregating sites = string length minus the leading g.
 newsites = NULL
 
 for (i in 1:nn) { #parse each g0101… line into a 0/1 vector and prepend c(1,1, …) as sentinels.
  gntp.s = as.numeric(strsplit(as.character(inpop$V1[i]),"")[[1]][2:(loci+1)])
  newsites = rbind(newsites,c(1,1,gntp.s))
 }
 
 
 out.pos = c(-2,-1,sample(sleng,loci)) #positions vector: sentinels at −2 and −1, then random unique positions for the M input SNPs.
 
 return(list(sites=newsites,pos=out.pos,fit=rep(1,nn/2)))

}


arg = commandArgs(TRUE) #parse command-line parameters: pur.sel (per-derived deleterious effect on nonsyn sites) and replicate rr.
pur.sel = as.numeric(arg[1])
rr = arg[2]
#core constants — Ne (diploid), sequence length, mu per bp, and mu per copy (mu.t = mu.bp * sleng)
ne = 2000
sleng = 20000
mu.bp = 5e-7
mu.t = mu.bp*sleng

nonsyn = rep(0,20000) #nonsyn mask (1) across the 20 kb region; design yields 8 kb nonsyn total.
nonsyn[3001:9000] = 1
nonsyn[12001:18000] = 1
nonsyn[(1:6000)*3] = 0
#8k nonsynonymous, 4k synonymous, 8k intron

nuc.st = rep(3,20000) #nuc.st: 1 = nonsyn, 2 = syn (3rd codon pos), 3 = intron; two 6-kb coding blocks.
nuc.st[3001:9000] = rep(c(1,1,2),2000)
nuc.st[12001:18000] = rep(c(1,1,2),2000)

pop = gen.read.onegene(paste("Round3.5_Initial/Round3.5_inieq_rep_",rr,".txt",sep="")) #read the neutral starting population file for this replicate.

print (paste("Pop read finish",arg))

  write("Gen Mean_Nonsyn_Muts Tpi Tpi_Nonsyn Tpi_Syn Tpi_Intron Tw Tw_Nonsyn Tw_Syn Tw_Intron",paste("Round3.5_Initial/Round3.5_inieq_pursel_",pur.sel,"_rep_",rr,"_selhistory.txt",sep="")) #write header for the selection-history time series file.

for (gen in 1:40000) { #main loop to 40,000 generations (mutation–selection balance warmup).
 temp.aa = c(FALSE,FALSE,as.logical(nonsyn[pop$pos[c(-2,-1)]])) #compute diploid deleterious burden at nonsyn sites and calculate multiplicative fitness
 del.ct = apply(pop$sites[1:ne,temp.aa],1,sum)+apply(pop$sites[(ne+1):(2*ne),temp.aa],1,sum)
 pop$fit = pop$fit*((1-pur.sel)^del.ct)
 pop = gen.onegene.norec(pop,ne) #evolve one generation (parent sampling + mutation + cleanup).
 
 if ((gen <= 1000) | (gen %% 100 == 0)) {
  theta.pi = NULL
  theta.w = NULL
  
  temp.fr = apply(pop$sites,2,mean) #overall theta_pi and theta_w (normalized by 20 kb)
  theta.pi[1] = sum(2*temp.fr*(1-temp.fr))/20
  theta.w[1] = (length(pop$pos)-2)/(20*sum(1/1:(2*ne)))
  
  temp.bl = c(FALSE,FALSE,(nuc.st[pop$pos[c(-2,-1)]]==1)) #class-specific theta_pi and theta_w for nonsyn, normalized by 8 kb.
  temp.fr = apply(pop$sites[,temp.bl],2,mean)
  theta.pi[2] = sum(2*temp.fr*(1-temp.fr))/8
  theta.w[2] = sum(temp.bl)/(8*sum(1/1:(2*ne))4)
  
  temp.bl = c(FALSE,FALSE,(nuc.st[pop$pos[c(-2,-1)]]==2)) #class-specific theta_pi and theta_w for syn, normalized by 4 kb.
  temp.fr = apply(pop$sites[,temp.bl],2,mean)
  theta.pi[3] = sum(2*temp.fr*(1-temp.fr))/4
  theta.w[3] = sum(temp.bl)/(4*sum(1/1:(2*ne)))
  
  temp.bl = c(FALSE,FALSE,(nuc.st[pop$pos[c(-2,-1)]]==3)) #class-specific theta_pi and theta_w for intron, normalized by 8 kb.
  temp.fr = apply(pop$sites[,temp.bl],2,mean)
  theta.pi[4] = sum(2*temp.fr*(1-temp.fr))/8
  theta.w[4] = sum(temp.bl)/(8*sum(1/1:(2*ne)))
  
  write(c(gen,mean(del.ct),theta.pi,theta.w),paste("Round3.5_Initial/Round3.5_inieq_pursel_",pur.sel,"_rep_",rr,"_selhistory.txt",sep=""),ncolumns=10,append=TRUE) #append a stats row: Gen Mean_Nonsyn_Muts Tpi_all … to the selhistory file.
 }
 
 if (gen %% 10000 == 0) { #every 10k generations, write full population snapshot (…_g10k.txt, …_g10k_pos.txt).
  write.table(pop$sites,paste("Round3.5_Initial/Round3.5_inieq_pursel_",pur.sel,"_rep_",rr,"_g",(gen/1000),"k.txt",sep=""),row.names=FALSE,col.names=FALSE)
  write(pop$pos,paste("Round3.5_Initial/Round3.5_inieq_pursel_",pur.sel,"_rep_",rr,"_g",(gen/1000),"k_pos.txt",sep=""),ncolumns=1)
 }
 
}
 


