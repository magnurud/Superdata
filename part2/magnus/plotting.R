#!/usr/bin/env Rscript
# Script for plotting relevant data in project    

#---------------------------------------------------------------------------------------
# Printing figures to file
NOPRINT <- TRUE # Set to FALSE to print to files
#NOPRINT <- FALSE # Set to FALSE to print to files

pathToFig <- "../Latex/figures/"
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

printfig("taskc1", NOPRINT=NOPRINT)
plot(x1$nthreads, x1$time, xlab="threads", ylab="time")
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

plotTimeVsProc <- function(t1){
    n <- unique(t1$n)
    timelim <- c(min(t1$time), max(t1$time))

    sub <- t1[t1$n == n[1], ]
    plot(sub$nproc, sub$time, ylim=timelim, pch='+', xlab = "processes", ylab="time")

    for (i in 2:length(n)) {
        sub <- t1[t1$n == n[i], ]
        points(sub$nproc, sub$time, col=i, pch='+')
    }
    grid()
    legend(x = "topright", as.character(n), pch = rep('+', length(n)),
           col = 1:3, bg = "white")
}

printfig("taskbTimeProc1", NOPRINT=NOPRINT)
plotTimeVsProc(t1)
off(NOPRINT)
printfig("taskbTimeProc2", NOPRINT=NOPRINT)
plotTimeVsProc(t2)
off(NOPRINT)

plotTimeVsThreadsTimesProc <- function(t1){
    n <- unique(t1$n)
    timelim <- c(min(t1$time), max(t1$time))

    sub <- t1[t1$n == n[1], ]
    plot(sub$nthreads*sub$nproc, sub$time, ylim=timelim, pch='+', xlab = "threads * processes", ylab="time")

    for (i in 2:length(n)) {
        sub <- t1[t1$n == n[i], ]
        points(sub$nthreads*sub$nproc, sub$time, col=i, pch='+')
    }
    grid()
    legend(x = "topright", as.character(n), pch = rep('+', length(n)),
           col = 1:3, bg = "white")
}

printfig("taskbTimeProcTimesThreads", NOPRINT=NOPRINT)
plotTimeVsThreadsTimesProc(t3)
off(NOPRINT)




# Convergence plots #---------------------------------------------

conv <- read.table("./convergence.txt")
names(conv) <- names(x1)

# Plot with different types, + O , and so on, to distinguish between them


