library(ggplot2)
##args <- commandArgs(trailingOnly = TRUE)

avdz_theo <- function(m,l){
    u = sqrt(m**2-l**2)
    avdz = 2*l/m**2 + l**2/u*besselI(u,1)*( l*besselK(u,1)/u - besselJ(l,0)/besselJ(l,1)*besselK(u,0))
    return(avdz)
}


colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
colfunc<-colorRampPalette(c("red","royalblue"))


main_dir <- c('/scratch/madorisi/results/polymers')
resfiles <- file.path( main_dir, 'results_L_Lz*.dat')
data_files <- Sys.glob(resfiles)



pow2 <- function(x){ 2**x }

df_complete <- data.frame()
df_log <- data.frame()
for (f in data_files){
    print(f)
    df_complete <- rbind(df_complete,
                         read.table(f))
    df_log <- rbind(df_log,
                    read.table(f)[sapply(seq(0,13),pow),])
}


colnames(df_log) <- c('L','Lz','avjmp')
colors <- rev( colfunc(length( unique(df_log$avjmp) )) )
L_vs_Lz <- ggplot(df_log) +
           geom_point(aes(x=L,y=Lz,colour=factor(avjmp)),size=1.) +
           geom_line(aes(x=L,y=Lz,colour=factor(avjmp)),size=.3) +
           scale_y_log10() + scale_x_log10() +
           coord_fixed(ratio = 1) +
           scale_colour_manual('average jump',values = colors) +
           theme_bw() +
           theme(legend.key=element_blank())
           
           
##gyration radius and end to end distance

df_ree_rg <- read.table('/scratch/madorisi/results/polymers/results_Ree_Rg.dat')
colnames(df_ree_rg) <- c('Rg2','Ree','totL','avjmp')


Rg_Ree_plot <- ggplot(df_ree_rg) +
           geom_point(aes(x=totL,y=sqrt(Rg2),colour=factor(avjmp)),shape=15,size=2.) +
           geom_point(aes(x=totL,y=Ree,colour=factor(avjmp)),shape=16,size=2.) +
           ##geom_line(aes(x=L,y=Lz,colour=factor(avjmp)),size=.3) +
           scale_y_log10() + scale_x_log10() +
           xlab('L') + ylab('Ree, Rg') +
           ##coord_fixed(ratio = 1) +
           scale_colour_manual('average jump',values = colors) +
           theme_bw() +
           theme(legend.key=element_blank())


df1 <-read.table('/scratch/madorisi/results/polymers/results_L_Lz_1.950.dat',header=TRUE)
df2 <-read.table('/scratch/madorisi/results/polymers/results_L_Lz_0.090.dat',header=TRUE)
avgjmp = 1.950
avg_dz_th <- avdz_theo(2./avgjmp,1.)##0.5*2.4048*avgjmp**2
histo_dz <- ggplot(df1) +
            geom_histogram(aes(avg_dz),bins=100) +
            geom_vline(aes(xintercept=mean(avg_dz)),color="blue", linetype="dashed", size=1) +
            geom_vline(aes(xintercept=avg_dz_th),color="red", size=1)

histo_dl <- ggplot(df1,aes(avg_dl)) +
            geom_histogram(aes(avg_dl),bins=100) +
            geom_vline(aes(xintercept=mean(avg_dl)),color="blue", linetype="dashed", size=1) +
            geom_histogram(data=df2,aes(avg_dl),bins=100)
            ##geom_vline(aes(xintercept=mean(avg_dl)),color="blue", linetype="dashed", size=1) +
            ##geom_vline(aes(xintercept=avgjmp),color="red", size=1)
