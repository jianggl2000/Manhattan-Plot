# manhattan plot using base graphics
manhattan = function(dataframe, cols=colors()[c(26,32,48,73,54,62,44,84,145,107,115,128,135,139,143,257,371,435,448,454,461,493,132)], ymin=0, ymax="max", cex.x.axis=0.7, limitchromosomes=1:26, suggestiveline=1:30, genomewideline=-log10(5e-8), annotate=NULL, ...) {
    cat("This function works only with plotting on the screen or pdf file.\n")
    d=dataframe
    if (any(limitchromosomes)) 
		d=d[d$CHR %in% limitchromosomes, ] #remove SNPs whose CHR outside of 1:26
	d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
	cat(nrow(d), "SNPs are to be plotted\n")
	
	labchr = sort(unique(d$CHR))
	labchr = ifelse(labchr==1,paste(labchr,"  ",sep=""),ifelse(labchr<=22, labchr, ifelse(labchr==23, "X", ifelse(labchr==24, "Y", ifelse(labchr==25,"XY", ifelse(labchr==26, "MT",labchr))))))

    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) 
		stop("Make sure your data frame contains columns CHR, BP, and P")

	d$logp = -log10(d$P)
	
    d$pos=NA
    ticks=NULL
    lastbase=0
    cols <- rep(cols, length=length(unique(d$CHR)))
    if (ymax=="max") 
		ymax<-ceiling(max(d$logp))
    if (ymin=="min")
		ymin=floor(min(d$logp))
		
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        #ticks=floor(length(d$pos))/2+1
        ticks = max(d$pos)/2+1
    } else {
        for (i in 1:numchroms) {
			mychr = sort(unique(d$CHR))[i]
        	if (i==1) {
    			d[d$CHR==mychr, ]$pos = d[d$CHR==mychr, ]$BP
				ticks = max(d[d$CHR==mychr,]$pos)/2
    		} else {
    			lastbase=lastbase+tail(subset(d,CHR==(unique(d$CHR)[i-1]))$BP, 1)
    			d[d$CHR==mychr, ]$pos = d[d$CHR==mychr, ]$BP+lastbase
				ticks = c(ticks, lastbase+tail(d[d$CHR==mychr,]$BP,1)/2)
    		}
    	}
    }
    # reset = readline("Would you accept the default 9*5 layout(y/n):")
    # if(reset=="y"){
        # devs = dev.list()
        # if(!is.null(devs))
            # dev.off()
        dev.new(width=9, height=5)
    # }
    par(xpd=F)
    par(mar=c(4.1, 4.1, 2.1, 2.1))
    if (numchroms==1) {
        with(d, plot(pos, logp, ylim=c(ymin,ymax), pch=".", ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), cex=0.4,cex.lab=1,cex.axis=1,las=1, ...))
    }	else {
        with(d, plot(pos, logp, ylim=c(ymin,ymax), pch=".",cex=2, ylab=expression(-log[10](italic(p))), xlab="", xaxt="n", type="n", cex.lab=1,cex.axis=0.8,las=1, ...))
        if (exists("suggestiveline"))
            abline(h=suggestiveline, col=colors()[353], lty=3)
        if (genomewideline)
            abline(h=genomewideline, col=colors()[134], lty=5)
        #axis(1, at=ticks, lab=labchr, las=2, lwd=0, lwd.ticks=1, line=-0.3, cex.axis=cex.x.axis)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, pch=".",cex=2, col=cols[icol], ...))
            icol=icol+1
        }
    }
   
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="green3", ...)) 
    }
    par(xpd=TRUE)
    halfchr = ceiling(length(labchr)/2)
    legend(-max(d$pos)/10,-ymax/20, paste("Chr",labchr[1:halfchr], sep=""), fill=cols[1:halfchr], horiz=T, cex=0.8, border=F,bty="n", x.intersp=0.3)
    legend(-max(d$pos)/10,-ymax/10, paste("Chr",labchr[(halfchr+1):length(labchr)],sep=""), fill=cols[(halfchr+1):length(labchr)], horiz=T, cex=0.8, border=F,bty="n", x.intersp=0.3)
    par(xpd=F)
}
