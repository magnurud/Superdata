#!/usr/bin/env Rscript
# Script for plotting relevant data in project    

#---------------------------------------------------------------------------------------
# Printing figures to file
#NOPRINT <- TRUE # Set to FALSE to print to files
NOPRINT <- FALSE # Set to FALSE to print to files

pathToFig <- "../Latex/Figures/"
printfig <- function(name, height = 6, pathFig = pathToFig, NOPRINT=FALSE){
    if (!NOPRINT){
        path <- paste(pathFig, name, '.pdf', sep = '')
        pdf(file = path, height = height)
        par(mai = c(0.6, 0.6, 0.1, 0.1)) # manipulate inner margins of subplots
        par(mgp = c(1.8, 0.7, 0)) # manipulate spacing of axis ticks, labels, text
        par(omi = c(0, 0, 0, 0)) # manipulate outer margin of full plot
    }
}
off <- function(NOPRINT=FALSE){
    if (!NOPRINT)
        invisible(dev.off())
}
#---------------------------------------------------------------------------------------


# Task C #--------------------------------------------
# One node
x1 <-  read.table("./taskc1.txt")
names(x1) <- c("nproc", "nthreads", "n", "time", "error")
x1NoMPI <- read.table("./taskc1NoMPI.txt")
names(x1NoMPI) <- names(x1)

printfig("taskc1", NOPRINT=NOPRINT)
plot(x1$nthreads, x1$time, xlab="threads", ylab="time")
points(x1NoMPI$nthreads[x1NoMPI$nthreads==12], x1NoMPI$time[x1NoMPI$nthreads==12], pch = 3, col = 2)
off(NOPRINT)


# Three nodes
x3 <-  read.table("./taskc3.txt")
names(x3) <- names(x1)
printfig("taskc2", NOPRINT=NOPRINT)
plot(x3$nthreads, x3$time, xlab="threads", ylab="time")
off(NOPRINT)




# Task B #-------------------------------------------
t1 <- read.table("./taskbTIME.txt")
t2 <- read.table("./taskbTIME2.txt")
t3 <- read.table("./taskbTIME3.txt")
names(t1) <- names(t2) <- names(t3) <- names(x1) 

plotSpeedupVsProc <- function(t1, ylim = NA, type = "speedup", vLines = NA,
                              legendLoc = "topleft", something = FALSE){
    n <- unique(t1$n)
    if (is.na(ylim[1])) {
        if (type == "efficiency") {
            ylim <- c(0, 1.2)
        }
        else if (type == "time")
            ylim <- c(min(t1$time), max(t1$time))
        else
            ylim <- c(0, max(t1$nproc))
    }

    if (something)
        t1$nproc <- t1$nproc*t1$nthreads
    sub <- t1[t1$n == n[1], ]
    sub <- sub[order(sub$nproc),]
    if (type == "efficiency") {
        if (something)
            plot(sub$nproc, sub$time[1]/sub$time/sub$nproc*2, ylim=ylim, pch='+', xlab = "processors", ylab="efficiency")
        else
            plot(sub$nproc, sub$time[1]/sub$time/sub$nproc, ylim=ylim, pch='+', xlab = "processors", ylab="efficiency")
    }
    else if (type == "time")
        plot(sub$nproc, sub$time, ylim=ylim, pch='+', xlab = "processors", ylab="time")
    else {
        if (something) 
            plot(sub$nproc, sub$time[1]/sub$time*2, ylim=ylim, pch='+', xlab = "processors", ylab="speedup")
        else
            plot(sub$nproc, sub$time[1]/sub$time, ylim=ylim, pch='+', xlab = "processors", ylab="speedup")
    }

    for (i in 2:length(n)) {
        sub <- t1[t1$n == n[i], ]
        sub <- sub[order(sub$nproc),]
        if (type == "efficiency") {
            if (something)
                points(sub$nproc, sub$time[1]/sub$time/sub$nproc*2, col=i, pch='+')
            else
                points(sub$nproc, sub$time[1]/sub$time/sub$nproc, col=i, pch='+')
        }
        else if (type == "time")
            points(sub$nproc, sub$time, col=i, pch='+')
        else {
            if (something) 
                points(sub$nproc, sub$time[1]/sub$time*2, col=i, pch='+')
            else
                points(sub$nproc, sub$time[1]/sub$time, col=i, pch='+')
        }
            
    }
    if (type == "efficiency") {
        abline(1, 0, lty=2)
    }
    else if (type == "speedup") {
        #if (something)
            #abline(0, 0.5, lty=2)
        #else
            abline(0, 1, lty=2)
    }
    else
        grid()

    if (!is.na(vLines[1])) {
        for (i in vLines) 
            abline(v=i, lty=3)
    }
    legend(x = legendLoc, as.character(n), pch = rep('+', length(n)),
           col = 1:3, bg = "white")
}

addLinearNetworkModel <- function(t1) {
    n   <- unique(t1$n)
    sub <- t1[t1$n == n[3], ]
    sub <- sub[order(sub$nproc),]

    tauC  <- 1e-6
    gamma <- 5e-9
    tauS  <- tauC + gamma*8*n[3]^2/sub$nproc
    Tlin  <- sub$time[3]/sub$nproc + tauS

    lines(sub$nproc, sub$time[3]/Tlin, lty = 2, col = 3)
}



printfig("taskbSpeedupProc1", NOPRINT=NOPRINT)
plotSpeedupVsProc(t1, ylim = c(0, 30), vLines=c(12, 24))
addLinearNetworkModel(t1)
off(NOPRINT)
printfig("taskbSpeedupProc2", NOPRINT=NOPRINT)
plotSpeedupVsProc(t2, ylim = c(0, 30), vLines=c(12, 24), something = TRUE)
off(NOPRINT)

plotSpeedupVsThreadsTimesProc <- function(t1, ylim = NA, type = "speedup", 
                                          vLines = NA, legendLoc = "topleft"){
    n <- unique(t1$n)
    if (is.na(ylim[1])) {
        if (type == "time")
            ylim <- c(min(t1$time), max(t1$time))
        else
            ylim <- c(0, max(t1$nproc*t1$nthreads))
    }

    sub <- t1[t1$n == n[1], ]
    sub <- sub[order(sub$nproc, sub$nthreads),]
    if (type == "time")
        plot(sub$nthreads*sub$nproc, sub$time, ylim=ylim, pch='+', xlab = "processors", ylab="time")
    else
        plot(sub$nthreads*sub$nproc, sub$time[1]/sub$time, ylim=ylim, pch='+', xlab = "processors", ylab="speedup")

    for (i in 2:length(n)) {
        sub <- t1[t1$n == n[i], ]
        sub <- sub[order(sub$nproc, sub$nthreads),]
        if (type == "time")
            points(sub$nthreads*sub$nproc, sub$time, col=i, pch='+')
        else
            points(sub$nthreads*sub$nproc, sub$time[1]/sub$time, col=i, pch='+')
    }
    if (type == "time")
        grid()
    else
        abline(0, 1, lty=2)

    if (!is.na(vLines[1])) {
        for (i in vLines) 
            abline(v=i, lty=3)
    }
    legend(x = legendLoc, as.character(n), pch = rep('+', length(n)),
           col = 1:3, bg = "white")
}

printfig("taskbSpeedupNodesTimesThreads", NOPRINT=NOPRINT)
#plotSpeedupVsThreadsTimesProc(t3, ylim = c(0, 25), vLines = c(12, 24))
plotSpeedupVsThreadsTimesProc(t3, ylim = c(0, 25))
off(NOPRINT)

#-------------------------
# Plotting speedup/(nproc) for 1 thread
printfig("taskbEfficiencyProc1", NOPRINT=NOPRINT)
plotSpeedupVsProc(t1, type = "efficiency", vLines=c(12, 24), legendLoc = "bottomleft")
off(NOPRINT)
printfig("taskbEfficiencyProc2", NOPRINT=NOPRINT)
plotSpeedupVsProc(t2, type = "efficiency", vLines=c(12, 24), legendLoc = "bottomleft", something = TRUE)
off(NOPRINT)

#------------------------
# Timing 

printfig("taskbTimeProc1", NOPRINT=NOPRINT)
plotSpeedupVsProc(t1, type = "time", legendLoc = "topright", vLines = c(12, 24))
off(NOPRINT)
printfig("taskbTimeProc2", NOPRINT=NOPRINT)
plotSpeedupVsProc(t2, type = "time", legendLoc = "topright", vLines = c(12, 24), something = TRUE)
off(NOPRINT)
printfig("taskbTimeNodesTimesThreads", NOPRINT=NOPRINT)
plotSpeedupVsThreadsTimesProc(t3, type = "time", legendLoc = "topright")
off(NOPRINT)


# Convergence plots #---------------------------------------------

conv <- read.table("./convergence.txt")
names(conv) <- names(x1)

# Plot with different types, + O , and so on, to distinguish between them
nproc <- unique(conv$nproc)
sub <- conv[conv$nproc == nproc[1],]
printfig("errVsn", NOPRINT=NOPRINT)
plot(sub$n, sub$error, pch = 1, col = 1, log="xy", xlab="n", ylab="error")
for (i in 2:length(nproc)) {
    sub <- conv[conv$nproc == nproc[i],]
    points(sub$n, sub$error, pch = i, col = i)
}
# Add line of slope -2
abline(max(conv$error), -2)

legend(x = "topright", as.character(nproc), pch = 1:length(nproc),
           col = 1:length(nproc), bg = "white")
off(NOPRINT)



# The rest of task d #-----------------------------------------------

# Timing as function of n^2
sub <- conv[conv$nproc == nproc[1],]
ylim <- c(min(conv$time/conv$n^2), max(conv$time/conv$n^2))
printfig("timeOverN2Vsn", NOPRINT=NOPRINT)
plot(sub$n, sub$time/sub$n^2, pch = 1, col = 1, xlab="n", 
     ylab="time/n^2", ylim = ylim)
for (i in 2:length(nproc)) {
    sub <- conv[conv$nproc == nproc[i],]
    points(sub$n, sub$time/sub$n^2, pch = i, col = i)
}
grid()
legend(x = "topright", as.character(nproc), pch = 1:length(nproc),
           col = 1:length(nproc), bg = "white")
off(NOPRINT)

# Timing as function of n^2
sub <- conv[conv$nproc == nproc[1],]
ylim <- c(min(conv$time/conv$n^2/log(conv$n)), 
          max(conv$time/conv$n^2/log(conv$n)))
printfig("timeOverN2LogNVsn", NOPRINT=NOPRINT)
plot(sub$n, sub$time/sub$n^2/log(sub$n), pch = 1, col = 1, xlab="n", 
     ylab="time/(n^2 log(n))", ylim = ylim)
for (i in 2:length(nproc)) {
    sub <- conv[conv$nproc == nproc[i],]
    points(sub$n, sub$time/sub$n^2/log(sub$n), pch = i, col = i)
}
grid()
legend(x = "topright", as.character(nproc), pch = 1:length(nproc),
           col = 1:length(nproc), bg = "white")
off(NOPRINT)
