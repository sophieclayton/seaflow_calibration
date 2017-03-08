library(flowCore)
library(splancs)
library(plotrix)
library(caroline)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))


plot.cytogram <- function (evtopp, para.x = "FSC.small.stuff", para.y = "X692.40.small.stuff", ...)
{
    cols <- colorRampPalette(c("blue4", "royalblue4", "deepskyblue3",
        "seagreen3", "yellow", "orangered2", "darkred"))
    par(pty = "s")
    plot(evtopp[, c(para.x, para.y)], pch = 16, cex = 0.3,col = densCols(log10(evtopp[, c(para.x, para.y)]),
          colramp = cols), xlim = c(1, 10^4), ylim = c(1, 10^4), log='xy',...)
          }

plot.vct.cytogram <- function (opp, para.x = "fsc_small", para.y = "chl_small", ...)
          {
        plot(opp[, c(para.x, para.y)], pch = 16, cex = 0.3, col = as.numeric(as.factor(opp$pop)),
                    xlim = c(1, 10^4), ylim = c(1, 10^4), log='xy', ...)
              legend("topleft", legend = (unique(opp$pop)), col = unique(as.numeric(as.factor(opp$pop))),
                  pch = 16, pt.cex = 0.6, bty = "n")
                    }


#################
## BATCH FILES ##
#################

# change this to point to the path where the .fcs files are
path <- paste("/Volumes/ceg/Sophie/bead_calibration/influx_bead_data")
savepath <- path
file.list <- dir(path, pattern = ".fcs$", recursive=F, full.names=T)

summary.table <- NULL

draw.gate <- TRUE

for (file in file.list) {
    print(paste("processing file:",file))


#file <- file.list[2]
###############
## read FCS ###
###############

fcs <- read.FCS(file, transformation=T, emptyValue=F)
opp <- tab2df(exprs(fcs))
opp$pop <- 0


##############
### GATING ###
##############

### NOISE & BEADS
x <- subset(opp, pop==0)
if(draw.gate) plot.cytogram(x, "X692.40.small.stuff", "X580.30", main="NOISE & BEADS & SYN")


print("Gating Beads")
if(draw.gate) poly.beads <- getpoly(quiet=TRUE)
beads <- subset(x,inout(x[,c("X692.40.small.stuff","X580.30")],poly=poly.beads, bound=TRUE, quiet=TRUE))
opp[row.names(beads),'pop'] <- "beads"

# print("Gating Synecho")
# if(draw.gate) poly.syn <- getpoly(quiet=TRUE)
# syn <- subset(x,inout(x[,c("X692.40.small.stuff","X580.30")],poly=poly.syn, bound=TRUE, quiet=TRUE))
# opp[row.names(syn),'pop'] <- "synecho"
#


### PRO & PicoEUK
x <- subset(opp, pop==0)
if(draw.gate) {
    plot.cytogram(x, "FSC.small.stuff", "X692.40.small.stuff", main="PRO & PicoEuk")
        points(beads$FSC.small.stuff, beads$X692.40.small.stuff, col="grey", pch=16, cex=0.4)
        }

# print("gating PRO")
# if(draw.gate) poly.pro <- getpoly(quiet=TRUE)
# pro <- subset(x,inout(x[,c("FSC.small.stuff","X692.40.small.stuff")],poly=poly.pro, bound=TRUE, quiet=TRUE))
# opp[row.names(pro),'pop'] <- "prochloro"


print("gating PICO")
if(draw.gate) poly.pico <- getpoly(quiet=TRUE)
pico <- subset(x,inout(x[,c("FSC.small.stuff","X692.40.small.stuff")],poly=poly.pico, bound=TRUE, quiet=TRUE))
opp[row.names(pico),'pop'] <- "picoeuk"

###################
### SAVE PLOT ###
###################
png(paste0(file,".png"),width=9, height=12, unit='in', res=100)

par(mfrow=c(2,2))
plot.vct.cytogram(opp, "FSC.small.stuff","X692.40.small.stuff")
plot.vct.cytogram(opp, "FSC.small.stuff","X580.30")
plot.vct.cytogram(opp, "X692.40.small.stuff","X580.30")
plot.vct.cytogram(opp, "SSC","X692.40.small.stuff")

dev.off()

###############
### SUMMARY ###
###############

stat.table <- NULL
for(i in unique(opp$pop)){
#print(i)
if(i == 0) next
p <- subset(opp, pop == i)
n <- nrow(p)
if(n ==0) {
fsc <- 0
chl <- 0
    }else{
fsc <- round(median(p$FSC.small.stuff))
chl <- round(median(p$X692.40.small.stuff))
ssc <- round(median(p$SSC))
pe <- round(median(p$X580.30))
}
var <- cbind(i,n,fsc,chl,ssc,pe)
stat.table <- rbind(stat.table, var)
}


table <- data.frame(cbind(stat.table, file=basename(file)))
summary.table <- rbind(summary.table, table)

}

write.csv(summary.table,file=paste(savepath,"/summary.csv", sep=""), row.names=FALSE)




path <- paste("~/Documents/DATA/SeaFlow/Calibration/FSCcalibration/Feb-23-2017")
path <- paste("~/Documents/DATA/SeaFlow/Calibration/FSCcalibration/Large_picoeuks")

savepath <- path
summary.table <- read.csv(file=paste(savepath,"/summary.csv", sep=""))
### FEB 10
#volume <- c(rep(50,4), rep(50,4), rep(50,4), rep(50,4),rep(150,4),rep(50,4),rep(100,8),rep(100,4),rep(50,4))
#name <- c("TAPS-1335", "E. hux","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus","Minutocellus","Navicula","TW")
### FEB 15
# volume <- c(rep(50,4), rep(100,4), rep(100,4), rep(50,4),rep(50,4),rep(100,4),rep(100,4),rep(100,4),rep(50,4),rep(50,4))
# name <- c("TAPS-1335", "E. hux","TW","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus","Minutocellus","Navicula")
# ### FEB 17
# volume <- c(rep(100,2), rep(150,2), rep(100,2), rep(100,2),rep(100,2),rep(200,2),rep(100,2),rep(100,2),rep(100,2))
# name <- c("TAPS-1335", "E. hux","TW","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus","Navicula")
# ### FEB 21
# volume <- c(rep(50,2), rep(150,2), rep(50,2), rep(50,2),rep(50,2),rep(150,2),rep(100,2),rep(100,2),rep(50,2))
# name <- c("TAPS-1335", "E. hux","TW","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus","Navicula")
# ### FEB 23
# volume <- c(rep(50,2), rep(100,2), rep(50,2), rep(50,2),rep(100,2),rep(50,2),rep(50,2),rep(50,2),rep(50,2),rep(50,2),rep(50,2))
# name <- c("TAPS-1335", "E. hux","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus-D","Minutocellus","Navicula-D","Navicula","TW")
# ### FEB 23 - POC
# volume <- c(rep(300.45,2), rep(300,2), rep(200.18,2), rep(200.1,2),rep(300.28,2),rep(302,2),rep(200.57,2),rep(199.97,2),rep(300.25,2),rep(300.24,2))
# name <- rep(c("Navicula","PT-632","TAPS-1335","TAPS-3367","TW-3365"),each=2)
# ### FEB 24
# volume <- rep(100,6)
# name <- c("E. hux","Licmophora","Micromonas")

### need to read in the volume run for each sample here

# summary.table$volume <- volume
# summary.table$abundance <- round(1000*summary.table$n / summary.table$volume)
# summary.table <- subset(summary.table, n > 0)




# beads <- subset(summary.table, i == 'beads')
# culture <- subset(summary.table, i == 'picoeuk')

# par(mar=c(4,20,2,2))
# barplot(culture$fsc/beads$fsc,horiz=T,names.arg=culture$file,las=1)
# mtext("Normalized Light scattering to 2 µm beads", 1, line =3)


# culture$abundance
