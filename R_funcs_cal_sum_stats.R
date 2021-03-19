sumstat.cal <- function(geno.both){
  icount <- 1
  sumstat = matrix(ncol = 95,nrow = 1)
  gt = make.single(geno.both)
  sumstat[icount,1] =  pairwise(geno.both)
  sumstat[icount,2] =  pairwise(geno.both[kola.chroms,])
  sumstat[icount,3] =  pairwise(geno.both[yamal.chroms,])
  sumstat[icount,4] =  pairwise(geno.both[kolyma.chroms,])
  sumstat[icount,5] =  pairwise(geno.both[guev.chroms,])
  sumstat[icount,6] =  sumstat[icount,2] / sumstat[icount,1]
  sumstat[icount,7] =  sumstat[icount,3] / sumstat[icount,1]
  sumstat[icount,8] =  sumstat[icount,4] / sumstat[icount,1]
  sumstat[icount,9] =  sumstat[icount,5] / sumstat[icount,1]
  sumstat[icount,10] = length(geno.both[1,])
  sumstat[icount,11] = sum(apply(geno.both[kola.chroms,],2,var) > 0)
  sumstat[icount,12] = sum(apply(geno.both[yamal.chroms,],2,var) > 0)
  sumstat[icount,13] = sum(apply(geno.both[kolyma.chroms,],2,var) > 0)
  sumstat[icount,14] = sum(apply(geno.both[guev.chroms,],2,var) > 0)
  sumstat[icount,15:25] = pairwise.sing(gt)
  sumstat[icount,26:36] = pairwise.sing(gt[kola,])
  sumstat[icount,37:47] = pairwise.sing(gt[yamal,])
  sumstat[icount,48:58] = pairwise.sing(gt[kolyma,])
  sumstat[icount,59:69] = pairwise.sing(gt[guev,])
  sumstat[icount,70:75] = summarise_fspec(fspec(gt))
  sumstat[icount,76:80] = cal.class(gt[kola,])
  sumstat[icount,81:85] = cal.class(gt[yamal,])
  sumstat[icount,86:90] = cal.class(gt[kolyma,])
  sumstat[icount,91:95] = cal.class(gt[guev,])
  return(sumstat)
}
cal.class  <- function(loci){
  nd <- nrow(loci) # nsamp  
  nnm1 <- nd / (nd - 1)
  
  p1 <- apply(loci, 2, sum) / nd
  segsite <- length(which(p1 != 0.0))
  pi <- sum(2 * p1 * (1 - p1) * nnm1)
  hfay<- (-1) * sum(2 * p1 * (2 * p1 - 1) / nnm1) 
  thetah <- 2 * sum( p1 ^ 2 ) / (nd * (nd - 1))
  # Taj'D 
  a1 <- sum(1 / seq(1:(nd-1)))
  a2 <- sum(1 / seq(1:(nd-1)) ^ 2)
  b1 <- (nd + 1) / (3 * (nd - 1))
  b2 <- (2* (nd * nd + nd + 3)) / (9 * nd * (nd - 1))
  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((nd + 2)/(a1 * nd)) + (a2 / (a1 ^ 2))
  e1 <- c1 / a1 
  e2 <- c2 / (a1 ^2 + a2 )
  tajd <- (pi - (segsite / a1)) / sqrt((e1 * segsite) + (e2 * segsite * (segsite - 1)))
  return(c(pi,segsite, tajd, thetah, hfay))
}

pairwise = function(m1)
{
	d1 = dist(m1,method="man")
	mean(d1) 
}

make.single = function(a)
{
	nhap = length(a[,1])
	nind = floor(nhap/2)
	if(nhap != 2*nind)stop("make.single: number of haplotypes is not even")
	nsnp = length(a[1,])
	h1 = a[seq(1,nhap,by=2),]
	h2 = a[seq(2,nhap,by=2),]
	s1 = h1 + h2
	s1
}

pairwise.sing = function(s)
{
	d1 = dist(s) 
	a1 = mean(d1)
	a2 = sd(d1)
	a3 = quantile(d1,c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99))
	as.numeric(c(a1,a2,a3))
	
}

fspec = function(s,wrap=T)
{
	nsnp = length(s[1,])
	nhap = length(s[,1])*2
	freqs = apply(s,2,sum)
	t1 = table(freqs)
	#return(t1)
	f1 = as.numeric(names(t1))
	c1 = as.numeric(t1)
	if(min(f1) == 0 || max(f1) == nhap){ 
		if(min(f1) == 0){
			f1 = f1[-1]
			c1 = c1[-1]	
		}
		if(max(f1) == nhap){
			f1 = f1[-length(f1)]
			c1 = c1[-length(c1)]	
		}
		if(length(f1)!=length(c1))stop("fspec: f1,c1 have different lengths")
		if(length(f1) > nhap-1)stop("fspec: f1,c1 too long")
	}
	f2 = c(1:(nhap-1))
	c2 = rep(0,nhap-1)
	c2[f1] = c1
	if(wrap){
		c3 = c2 + rev(c2)
		c3[nhap/2] = c3[nhap/2]/2
		return(c3[1:(nhap/2)])
	}
	return(c2)
	
}

summarise_fspec = function(v)
{
	ip = length(v)
	if(ip <= 10)stop("summarise_fspec: currently assumes spectrum has length greater than 10")
	q1 = ceiling(rep(ip,6)/c(ip,10,5,2.5,1.75,1.25))
	a1 = cumsum(v)/sum(v)
	a1[q1]
}


getcovs = function(m1)
{

	xx = cov(m1)
	utri = upper.tri(xx)
	c1 = abs(xx[utri])
	c(mean(c1),sd(c1))
}


