setClass(Class = "sample",
	representation(	sample.sizes = "numeric",
					number.of.loci = "numeric",
					nbr.of.sampled.demes = "numeric",
					counts = "matrix")
)

setClass(Class = "likelihood",
	representation(	max.pi = "numeric",
					max.M = "numeric",
					pi.limits = "numeric",
					M.limits = "numeric")
)

SimulateCoalescentTree <- function(total.number.of.demes,effective.size,migration.rate,mutation.rate,number.of.sampled.demes,sample.sizes,number.of.loci) {
	.C('SimulateCoalescentTree',
	as.integer(total.number.of.demes),
	as.integer(effective.size),
	as.double(migration.rate),
	as.double(mutation.rate),
	as.integer(number.of.sampled.demes),
	as.integer(sample.sizes),
	as.integer(number.of.loci),
	PACKAGE = 'PopGenInfeR')
}

SimulateIslandModel <- function(number.of.simulations,mutation.model,total.number.of.demes,number.of.loci,number.of.sampled.demes,sample.sizes) {
	.C('SimulateIslandModel',
	as.integer(number.of.simulations),
	as.integer(mutation.model),
	as.integer(total.number.of.demes),
	as.integer(number.of.loci),
	as.integer(number.of.sampled.demes),
	as.integer(sample.sizes),
	PACKAGE = 'PopGenInfeR')
}

migraine.colors <- function (n = 64,redshift = 1) { ## added redshift to tim.colors to control the amount of red...
	orig <- c(	"#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
				"#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
				"#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
				"#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
				"#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
				"#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
				"#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
				"#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
				"#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
				"#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
				"#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
				"#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
				"#AF0000", "#9F0000", "#8F0000", "#800000")
	orig[1:20] <- topo.colors(64)[1:20] ## the blues in tim.colors are too dark for filled.contours... ?pie to test other combinations
	if (n == 64 && redshift == 1) return(orig)
	rgb.tim <- t(col2rgb(orig))
	temp <- matrix(NA, ncol = 3, nrow = n)
	x <- (seq(0, 1, , 64)^(redshift)) ## values in [0,1]
	xg <- seq(0, 1, , n)
	for (k in 1:3) {
		hold <- fields::splint(x, rgb.tim[, k], xg)
		hold[hold < 0] <- 0
		hold[hold > 255] <- 255
		temp[, k] <- round(hold)
	}
	rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}

log.beta.binomial <- function(n,k,M,pi) {
	result <- lgamma(M) - lgamma(M + n) + lchoose(n,k) + lgamma(M * pi + k) - lgamma(M * pi) + lgamma(M * (1 - pi) + (n - k)) - lgamma(M * (1 - pi))
	return(result)
}

sim.inference.model <- function(number.of.sampled.demes,sample.sizes,M,pi) {
	number.of.loci <- length(pi)
	counts <- matrix(nrow = number.of.loci,ncol = number.of.sampled.demes)
	X <- vector(mode = "numeric",length = number.of.sampled.demes)

	for (i in 1:number.of.loci) {
		X <- rbeta(number.of.sampled.demes,M * pi[i],M * (1 - pi[i]))
		for (j in 1:number.of.sampled.demes) {
			counts[i,j] <- rbinom(1,sample.sizes,X[j])
		}
	}
	sample <- new("sample")
	sample@number.of.loci <- number.of.loci
	sample@nbr.of.sampled.demes <- number.of.sampled.demes
	sample@sample.sizes <- sample.sizes
	sample@counts <- counts
	return(sample)		
}

expected.fst <- function(total.number.of.demes,effective.size,migration.rate,mutation.rate) {
	k <- 2
	mu <- mutation.rate
	N <- effective.size
	m <- migration.rate
	n <- total.number.of.demes
	gamma <- (1 - mu * k / (k - 1))^2
	a <- (1 - m)^2 + m^2 / (n - 1)
	b <- (1 - a) / (n - 1)
	d <- (a - b)
	result <- (gamma * d) / (gamma * d + N * (1 - gamma * d)) 
	return(result)
}

fst <- function(sample) {
	n <- sample@nbr.of.sampled.demes * sample@sample.sizes
	nsqr <- sample@nbr.of.sampled.demes * sample@sample.sizes^2
	nc <- (n - nsqr / n) / (sample@nbr.of.sampled.demes - 1)
	counts <- sample@counts
	frequencies <- counts / sample@sample.sizes
	marginal.frequencies <- rowSums(counts) / n
	SSI <- rowSums((frequencies - frequencies^2) * sample@sample.sizes) + rowSums(((1 - frequencies) - (1 - frequencies)^2) * sample@sample.sizes)
	SSP <- rowSums((frequencies - marginal.frequencies)^2 * sample@sample.sizes) + rowSums(((1 - frequencies) - (1 - marginal.frequencies))^2 * sample@sample.sizes)
	MSI <- SSI / (n - sample@nbr.of.sampled.demes) 
	MSP <- SSP / (sample@nbr.of.sampled.demes - 1)
	result <- sum(MSP - MSI) / sum(MSP + (nc - 1) * MSI)
	return(result)	
}
 
sim.coalescent <- function(total.number.of.demes,effective.size,migration.rate,mutation.rate,number.of.sampled.demes,sample.sizes,number.of.loci) {
	SimulateCoalescentTree(total.number.of.demes,effective.size,migration.rate,mutation.rate,number.of.sampled.demes,sample.sizes,number.of.loci)	
	counts <- as.matrix(read.table("tmp.dat"))
	file.remove("tmp.dat")
	sample <- new("sample")
	sample@number.of.loci <- number.of.loci
	sample@nbr.of.sampled.demes <- number.of.sampled.demes
	sample@sample.sizes <- sample.sizes
	sample@counts <- counts
	return(sample)		
}

maximum.likelihood <- function(sample,alpha = 0.05,M = seq(0.1,9.99,0.1),pi = seq(0.01,0.99,0.01),graphics = TRUE,true.M,true.pi) {
	if (sample@number.of.loci > 1) {
		stop("Only the samples with a single locus can be analyzed with this function...")
	}
	log.likelihood <- function(M,pi) {
		n <- sample@sample.sizes
		result <- 0
		for (i in 1:sample@number.of.loci) {
			for (j in 1:sample@nbr.of.sampled.demes) {
				k <- sample@counts[i,j]
				result <- result + log.beta.binomial(n,k,M,pi)
			}
		}
		return(result)
	}
	z <- outer(M,pi,log.likelihood)																			# Matrix of log(likelihood) with M rows and pi columns
	profile.pi <- apply(z,2,max)																			# Likelihood profile for pi
	profile.M <- apply(z,1,max)																				# Likelihood profile for M
	likelihood.ratio <- exp(z - max(z))
	likelihood.ratio.pi <- exp(profile.pi - max(profile.pi))
	likelihood.ratio.M <- exp(profile.M - max(profile.M))
	likelihood.ratio.pvalue <- pchisq(-2 * (z - max(z)),df = 2,lower.tail = FALSE)
	likelihood.ratio.pi.pvalue <- pchisq(-2 * (profile.pi - max(profile.pi)),df = 1,lower.tail = FALSE)
	likelihood.ratio.M.pvalue <- pchisq(-2 * (profile.M - max(profile.M)),df = 1,lower.tail = FALSE)
	max.pi.idx <- which(profile.pi == max(profile.pi), arr.ind=TRUE)
	max.M.idx <- which(profile.M == max(profile.M), arr.ind=TRUE)
	pi.limits <- range(pi[likelihood.ratio.pi.pvalue >= alpha / 2])
	M.limits <- range(M[likelihood.ratio.M.pvalue >= alpha / 2])
	if (graphics) {
		dev.new()
		plot.new()
# Top left:
		par(new = "TRUE",plt = c(0.925,0.945,0.125,0.4),cex.axis = 0.6,las = 1,mgp = c(2.5,0.75,0))
		filled.legend(M,pi,z,color.palette = migraine.colors)
		par(new = "TRUE",plt = c(0.1,0.45,0.625,0.9),cex.axis = 0.83, cex.lab = 0.83,mgp = c(2.5,1,0))
		surface <- persp(M,pi,z,theta = 30,phi = 30,expand = 0.5,col = "lightblue",ltheta = 60,shade = 0.25,ticktype = "detailed",xlab = "M",ylab = "pi",zlab = "",main = "Log(likelihood)")
		points(trans3d(M[max.M.idx],pi[max.pi.idx],max(z),pmat = surface),col = "red",pch = 3)
		if (!missing(true.M) & !missing(true.pi)) {
			points(trans3d(true.M,true.pi,max(z),pmat = surface),col = "red",pch = 1)
		}
# Top right:
		par(new = "TRUE",plt = c(0.6,0.95,0.625,0.9),cex.axis = 0.83, cex.lab = 0.83,las = 0,mgp = c(2.5,1,0))
		plot(pi,profile.pi,xlab = expression(pi),ylab = "Log(likelihood)",type = "l",main = bquote(hat(pi) == .(pi[max.pi.idx])))
 		points(pi[max.pi.idx],profile.pi[max.pi.idx],col = "red",pch = 3)
		if (!missing(true.pi)) {
			points(true.pi,profile.pi[max.pi.idx],col = "red",pch = 1)
		}
 		abline(v = pi[max.pi.idx],lty = 3,col = "red")
# Bottom left:
		par(new = "TRUE",plt = c(0.1,0.45,0.125,0.4),cex.axis = 0.83, cex.lab = 0.83,las = 0,mgp = c(2.5,1,0))
		plot(M,profile.M,xlab = bquote(.(italic(M)~"= 2"~italic(N)[e]~italic(m))),ylab = "Log(likelihood)",type = "l",main = bquote(hat(italic(M)) == .(M[max.M.idx])))
 		points(M[max.M.idx],profile.M[max.M.idx],col = "red",pch = 3)
		if (!missing(true.M)) {
			points(true.M,profile.M[max.M.idx],col = "red",pch = 1)
		}
 		abline(v = M[max.M.idx],lty = 3,col = "red")
# Bottom right:
		par(new = "TRUE",plt = c(0.6,0.9,0.125,0.4),cex.axis = 0.83, cex.lab = 0.83,las = 0,mgp = c(2.5,1,0))
		filled.contour3(M,pi,z,color.palette = migraine.colors,main = "Log(likelihood)",xlab = bquote(.(italic(M)~"= 2"~italic(N)[e]~italic(m))),ylab = expression(pi))
		contour(M,pi,z,xlab = "M",ylab = "pi",col = "white",add = TRUE)
		points(M[which(z == max(z), arr.ind=TRUE)[1]],pi[which(z == max(z), arr.ind=TRUE)[2]],col = "white",pch = 3,lwd = 1)
		if (!missing(true.M) & !missing(true.pi)) {
			points(true.M,true.pi,col = "white",pch = 1)
		}
		dev.new()
		plot.new()
# Top left: 
		par(new = "TRUE",plt = c(0.425,0.445,0.625,0.9),cex.axis = 0.6,las = 1,mgp = c(2.5,0.75,0))
		filled.legend(M,pi,likelihood.ratio,color.palette = migraine.colors,zlim = c(0,1))
		par(new = "TRUE",plt = c(0.1,0.4,0.625,0.9),cex.axis = 0.83, cex.lab = 0.83,mgp = c(2.5,1,0))
		filled.contour3(M,pi,likelihood.ratio,color.palette = migraine.colors,main = "Likelihood ratio", xlab = bquote(.(italic(M)~"= 2"~italic(N)[e]~italic(m))),ylab = expression(pi))
		points(M[which(z == max(z), arr.ind=TRUE)[1]],pi[which(z == max(z), arr.ind=TRUE)[2]],col = "white",pch = 3,lwd = 1)
		contour(M,pi,likelihood.ratio,xlab = "M",ylab = "pi",col = "white",add = TRUE)
		if (!missing(true.M) & !missing(true.pi)) {
			points(true.M,true.pi,col = "white",pch = 1)
		}
# Top right: 
		par(new = "TRUE",plt = c(0.6,0.95,0.625,0.9),cex.axis = 0.83, cex.lab = 0.83,las = 0,mgp = c(2.5,1,0))
		plot(pi,likelihood.ratio.pi,xlab = expression(pi),ylab = "Likelihood ratio",type = "l",main = bquote(CI[.(1 - alpha)] == ~.("[")~.(pi.limits[1])~","~.(pi.limits[2])~"]"),log = "y")
  		points(pi[max.pi.idx],max(likelihood.ratio.pi),col = "red",pch = 3)
		if (!missing(true.pi)) {
			points(true.pi,max(likelihood.ratio.pi),col = "red",pch = 1)
		}
		abline(v = pi[max.pi.idx],lty = 3,col = "red")
 		abline(v = pi.limits,lty = 3,col = "green")
		abline(h = exp(-qchisq(p = alpha,df = 1,lower.tail = FALSE) / 2),lty = 3,col = "blue")
# Bottom left:
		par(new = "TRUE",plt = c(0.1,0.45,0.125,0.4),cex.axis = 0.83, cex.lab = 0.83,las = 0,mgp = c(2.5,1,0))
		plot(M,likelihood.ratio.M,xlab = bquote(.(italic(M)~"= 2"~italic(N)[e]~italic(m))),ylab = "Likelihood ratio",type = "l",main = bquote(CI[.(1 - alpha)] == ~.("[")~.(M.limits[1])~","~.(M.limits[2])~"]"),log = "y")
 		points(M[max.M.idx],max(likelihood.ratio.M),col = "red",pch = 3)
		if (!missing(true.M)) {
			points(true.M,max(likelihood.ratio.M),col = "red",pch = 1)
		}
 		abline(v = M[max.M.idx],lty = 3,col = "red")
 		abline(v = M.limits,lty = 3,col = "green")
 		abline(h = exp(-qchisq(p = alpha,df = 1,lower.tail = FALSE) / 2),lty = 3,col = "blue")
# Bottom right:
		par(new = "TRUE",plt = c(0.6,0.95,0.125,0.4),cex.axis = 0.83, cex.lab = 0.83,las = 0,mgp = c(2.5,1,0))
		filled.contour3(M,pi,likelihood.ratio.pvalue,color.palette = migraine.colors,main = "2D confidence interval", xlab = bquote(.(italic(M)~"= 2"~italic(N)[e]~italic(m))),ylab = expression(pi))
		points(M[which(z == max(z), arr.ind=TRUE)[1]],pi[which(z == max(z), arr.ind=TRUE)[2]],col = "white",pch = 3,lwd = 1)
		if (!missing(true.M) & !missing(true.pi)) {
			points(true.M,true.pi,col = "white",pch = 1)
		}
		contour(M,pi,likelihood.ratio.pvalue,xlab = "M",ylab = "pi",col = "white",levels = c(alpha),labels = c(paste(toString((1 - alpha) * 100),"%",sep = "")),add = TRUE)
	}	
	likelihood <- new("likelihood")
	likelihood@max.pi <- pi[max.pi.idx] 
	likelihood@max.M <- M[max.M.idx] 
	likelihood@pi.limits <- pi.limits
	likelihood@M.limits <- M.limits
	return(likelihood)
}

run.mcmc <- function(sample,chain.length = 1000,burnin = 0,range.M = c(0.001,10),delta.M = 2.5,delta.pi = 0.25,alpha = 0.05,graphics = TRUE,true.M,true.pi) {
	mcmc.chain <- matrix(ncol = (1 + sample@number.of.loci),nrow = chain.length)
	M <- runif(1) * (range.M[2] - range.M[1]) + range.M[1]
	pi <- vector(mode = "numeric",length = sample@number.of.loci)
	for (i in 1:sample@number.of.loci) {
		pi[i] <- runif(1)
	}
	LIK.1 <- lgamma(M) - lgamma(M + sample@sample.sizes)
	LIK.2 <- matrix(0,nrow = sample@number.of.loci,ncol = sample@nbr.of.sampled.demes)
	for (i in 1:sample@number.of.loci) {
		tmp <- lgamma(M * pi[i]) + lgamma(M * (1 - pi[i]))
		for (j in 1:sample@nbr.of.sampled.demes) {
			LIK.2[i,j] <- lgamma(M * pi[i] + sample@counts[i,j]) + lgamma(M * (1 - pi[i]) + sample@sample.sizes - sample@counts[i,j]) - tmp
		}
	}
	mcmc.chain[1,] <- c(M,pi)
	for (step in 2:chain.length) {
		new <- M + runif(1) * delta.M - delta.M / 2											# Update M
		if (new < range.M[1]) {
			new <- range.M[1] + (range.M[1] - new)
		}
		if (new > range.M[2]) {
			new <- range.M[2] - (new - range.M[2])				
		}
		prop.1 <- lgamma(new) - lgamma(new + sample@sample.sizes)
		prop.2 <- LIK.2 * 0
		for (i in 1:sample@number.of.loci) {
			tmp <- lgamma(new * pi[i]) + lgamma(new * (1 - pi[i]))
			for (j in 1:sample@nbr.of.sampled.demes) {
				prop.2[i,j] <- lgamma(new * pi[i] + sample@counts[i,j]) + lgamma(new * (1 - pi[i]) + sample@sample.sizes - sample@counts[i,j]) - tmp
			}
		}
        h <- sample@number.of.loci * sample@nbr.of.sampled.demes * (prop.1 - LIK.1)  + (sum(prop.2) - sum(LIK.2))            
		if (log(runif(1)) < h) {
			M <- new
			LIK.1 <- prop.1
			LIK.2 <- prop.2
		}
		mcmc.chain[step,1] <- M
    	for (i in 1:sample@number.of.loci) {
			new <- abs(pi[i] + runif(1) * delta.pi - delta.pi / 2)
  			if (new > 1) {
  				new <- 2 - new
  			}
			prop.2 <- LIK.2[i,] * 0 - lgamma(M * new) - lgamma(M * (1 - new))
			for (j in 1:sample@nbr.of.sampled.demes) {
				prop.2[j] <- prop.2[j] + lgamma(M * new + sample@counts[i,j]) + lgamma(M * (1 - new) + sample@sample.sizes - sample@counts[i,j])
			}
			h <- sum(prop.2) - sum(LIK.2[i,])										# Compute Metropolis ratio
			if (log(runif(1)) < h) {
				pi[i] <- new
				LIK.2[i,] <- prop.2
			}
			mcmc.chain[step,(i + 1)] <- pi[i]
		}
	}
	mcmc.chain <- mcmc.chain[-c(1:burnin),]
	if (graphics) {
		dev.new()
		if (sample@number.of.loci > 1) {
			tmp <- {"M"}
			for (i in 1:sample@number.of.loci) {
				 tmp <- append(tmp,paste("pi [locus ",toString(i),"]",sep = ""))
			}
			colnames(mcmc.chain) <- tmp
			plot(coda::as.mcmc(mcmc.chain))	
		} else {
			par(mfrow = c(3,2))
			z <- MASS::kde2d(mcmc.chain[,1],mcmc.chain[,2],lims = c(range.M,c(0,1)))
			surface <- persp(z,theta = 30,phi = 30,expand = 0.5,col = "lightblue",ltheta = 60,shade = 0.25,ticktype = "detailed",xlab = "M",ylab = "pi",zlab = "",main = "Posterior density")
			points(trans3d(z$x[which(z$z == max(z$z),arr.ind = TRUE)[1]],z$y[which(z$z == max(z$z),arr.ind = TRUE)[2]],max(z$z),pmat = surface),col = "red",pch = 3)
			
			if (!missing(true.M) & !missing(true.pi)) {
				points(trans3d(true.M,true.pi,max(z$z),pmat = surface),col = "red",pch = 1)
			}
			filled.contour3(z,color.palette = migraine.colors,main = "Posterior density",xlab = bquote(.(italic(M)~"= 2"~italic(N)[e]~italic(m))),ylab = expression(pi))
			contour(z,col = "white",levels = c(alpha),labels = c(paste(toString((1 - alpha) * 100),"%",sep = "")),add = TRUE)
			points(z$x[which(z$z == max(z$z),arr.ind = TRUE)[1]],z$y[which(z$z == max(z$z),arr.ind = TRUE)[2]],col = "white",pch = 3,lwd = 1)
			if (!missing(true.M) & !missing(true.pi)) {
				points(true.M,true.pi,col = "white",pch = 1)
			}
			pi.limits <-  coda::HPDinterval(coda::as.mcmc(mcmc.chain[,2]),prob = (1 - alpha))
			fit <- locfit::locfit(~mcmc.chain[,2])
			max.pi.idx <- which(fitted(fit) == max(fitted(fit)))
			plot(mcmc.chain[,2],cex = 0.075,ylim = c(0,1),xlab = "Iterations",ylab = expression(pi))
			plot(fit,get.data = TRUE,xlim = c(0,1),xlab = bquote(pi),ylab = "Posterior density",main = bquote(hat(pi) == .(mcmc.chain[max.pi.idx,2])~"   "~CI[.(1 - alpha)] == ~.("[")~.(pi.limits[1])~","~.(pi.limits[2])~"]"))
 			points(mcmc.chain[max.pi.idx,2],rep(max(fitted(fit)),length(max.pi.idx)),col = "red",pch = "+")
			if (!missing(true.pi)) {
				points(true.pi,max(fitted(fit)),col = "red",pch = 1)
			}
 			abline(v = mcmc.chain[max.pi.idx,2],lty = 3,col = "red")
 			abline(v = pi.limits,lty = 3,col = "green")
 			M.limits <- coda::HPDinterval(coda::as.mcmc(mcmc.chain[,1]),prob = (1 - alpha))
			fit <- locfit::locfit(~mcmc.chain[,1])
			max.M.idx <- which(fitted(fit) == max(fitted(fit)))
			plot(mcmc.chain[,1],cex = 0.075,ylim = range.M,xlab = "Iterations",ylab = bquote(.(italic(M)~"= 2"~italic(N)[e]~italic(m))))
			plot(fit,get.data = TRUE,xlim = range.M,xlab = bquote(.(italic(M)~"= 2"~italic(N)[e]~italic(m))),ylab = "Posterior density",main = bquote(hat(italic(M)) == .(mcmc.chain[max.M.idx,1])~"   "~CI[.(1 - alpha)] == ~.("[")~.(M.limits[1])~","~.(M.limits[2])~"]"))
 			points(mcmc.chain[max.M.idx,1],rep(max(fitted(fit)),length(max.M.idx)),col = "red",pch = "+")
 			if (!missing(true.M)) {
				points(true.M,max(fitted(fit)),col = "red",pch = 1)
			}
 			abline(v = mcmc.chain[max.M.idx,1],lty = 3,col = "red")
 			abline(v = M.limits,lty = 3,col = "green")
		}
	}  	
  	return(coda::as.mcmc(mcmc.chain))
}

generate.prior <- function(number.of.simulations,prior.theta,min.theta,max.theta,prior.M,min.M,max.M) {
	if (prior.theta == "UNI") {
		theta <- runif(number.of.simulations,min = min.theta,max = max.theta)
	}
	else if (prior.theta == "LOG") {
		theta <- exp(runif(number.of.simulations,min = log(min.theta),max = log(max.theta)))
	}
	else {
		stop("Unknown prior for theta (must be one of 'UNI' or 'LOG')...")
	}
	if (prior.M == "UNI") {
		M <- runif(number.of.simulations,min = min.M,max = max.M)
	}
	else if (prior.M == "LOG") {
		M <- exp(runif(number.of.simulations,min = log(min.M),max = log(max.M)))
	}
	else {
		stop("Unknown prior for theta (must be one of 'UNI' or 'LOG')...")
	}
	prior <- as.matrix(cbind(theta,M))
	write.table(file = "prior_val.in",prior,col.names = FALSE,row.names = FALSE,quote = FALSE)
	return(data.frame(prior))
}

sim.island.model <- function(number.of.simulations,mutation.model,total.number.of.demes,number.of.loci,number.of.sampled.demes,sample.sizes) {
	if (mutation.model == "IAM") {
		mutation.model = 0
	}
	else if (mutation.model == "KAM") {
		mutation.model = 1
	}
	else if (mutation.model == "SMM") {
		mutation.model = 2
	}
	else {
		stop("Unknown mutation model (must be one of 'IAM', 'KAM', or 'SMM')...")
	}
	SimulateIslandModel(number.of.simulations,mutation.model,total.number.of.demes,number.of.loci,number.of.sampled.demes,sample.sizes)
	summary_stats <- as.matrix(read.table("sim_stats.out",header = TRUE))
	file.remove("prior_val.in")
	file.remove("sim_stats.out")
	return(data.frame(summary_stats))
}

prcp.ca <- function(stats,target,...) {
	scaled.stats <- scale(stats)
	scaled.target <- (target - colMeans(stats)) / apply(stats,2,sd)
	pca <- ade4::dudi.pca(scaled.stats,scannf = FALSE,nf = 2)
	pca.target <- t(pca$c1) %*% t(scaled.target)
	percentage <- pca$eig / sum(pca$eig) * 100
	plot(pca$li,pty = 1,cex = 0.1,xlab = paste("Dimension 1 (",format(percentage[1],digits = 3),"%)",sep = ""),ylab = paste("Dimension 2 (",format(percentage[2],digits = 3),"%)",sep = ""),main = "Principal Compenent Analysis",cex.lab = 1.25,cex.main = 1.5,...)
	points(pca.target[1],pca.target[2],col = "red",pch = 16,cex = 1.75)
	return(pca)
}

