# This script includes the R code used in the analyses presented in:
#
# 	Zapata F, Villarroel D. A new species of Escallonia (Escalloniaceae) from
#	the inter-Andean tropical dry forests of Bolivia
#
#
# Code sections: preliminaries,  data organization, analyses, and figures.
# Each subsection is separated by -------------------------
#

-------------------------
#########################
#                       #
#    PRELIMINARIES      #
#                       #
#########################

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(labdsv)
library(vegan)
library(spdep)
library(ellipse)
library(mvtnorm)
library(dplyr)
library(ggmap)
library(maptools)
library(rgeos)
library(legendMap)


data = read.table("data.csv",
		   sep =",",
		   header = TRUE)

specimens = read.table("specimens.csv",
			sep = ",",
			header = TRUE,
			quote = "",
			fill = TRUE,
			row.names = NULL)


-------------------------
###########################
#                         #
#    DATA ORGANIZATION    #
#                         #
###########################

data = arrange(data, species)
species_names = as.character(unique(data$species))
species_names2 = c( species_names[11], species_names[1:10],  species_names[12:39])
species_names3 = species_names2
species_names3[1] = "New species"

harrisii = subset(data, species == "harrisii")

out = c( "harrisii",
         "callcottiae",
         "gayana",
         "ledifolia",
         "myrtilloides",
         "polifolia",
         "rosea",
         "serrata" )

out_names = out
out_names[1] = "New species"


# MAPS

southam = map_data(map = "world",
		   region = c( "argentina",
			      "bolivia",
			      "brazil",
			      "chile",
			      "colombia",
			      "ecuador",
			      "guyana",
			      "paraguay",
			      "peru",
			      "suriname",
			      "uruguay",
			      "venezuela",
			      "panama",
			      "costa rica" ))

bolivia = map_data(map = "world",
		   region = c( "bolivia" ))

subdata = data[data$species=="harrisii",]

zoom_box = data.frame(xmin = -66,
		      xmax = -62,
		      ymin = -22,
		      ymax = -18 )

eh_bbox = make_bbox( lon = subdata$long, lat = subdata$lat, f = 0.1 )

map2 = get_map( location = eh_bbox,
	       maptype = "satellite",
	       source = "google",
	       color = "bw",
	       zoom = 8 )



-------------------------
#######################
#                     #
#    DATA ANALYSES    #
#                     #
#######################

# These analyses use the scripts to detect discontinuities in morphology presented in:

# Zapata & Jiménez (2012). Species delimitation: Inferring gaps in morphology across geography.
# Syst. Biol. 61:179-194

# DISCONTINUITIES IN LEAF MORPHOLOGY

# Are there discontinuities in the bivariate space defined by lamina length vs.lamina width?

# Lists to store plots
manifold_plots = list()
mixture_plots = list()
proportion_plots = list()
leaf_plots = list()

#Gamma level of statistical confidence to calculate statistical tolerance regions (see Krishnamoorty and Mathew 1999).
gamma = 0.90

#Frequency cutoff
freq_cutoff = 0.85

#Spp in mixture
p = 2

#Compare all species vs. new species (E. harrisii)

for (i in unique(data$species)[unique(data$species) != c("harrisii")])
    {
        subdata = data[data$species %in% c(paste(i), "harrisii"),]
        SPA <- "harrisii"
        SPB <- paste(i)

        n_SPA <- length(subdata$LAMLEN[subdata$species==SPA])
        m1 <- array(c(mean(subdata$LAMLEN[subdata$species==SPA]), mean(subdata$LAMWID[subdata$species==SPA])), dim=c(2,1))
        Z1 <- var(cbind(subdata$LAMLEN[subdata$species==SPA], subdata$LAMWID[subdata$species==SPA]))

        n_SPB <- length(subdata$LAMLEN[subdata$species==SPB])
        m2 <- array(c(mean(subdata$LAMLEN[subdata$species==SPB]), mean(subdata$LAMWID[subdata$species ==SPB])), dim=c(2,1))
        Z2 <- var(cbind(subdata$LAMLEN[subdata$species==SPB], subdata$LAMWID[subdata$species==SPB]))

        #Estimate ridgeline manifold
        x <- c()
        y <- c()
        alfa <- seq(0,1,0.001)
        for (n in 1:length(alfa))
            {
                a <- solve((1-alfa[n])*solve(Z1) + alfa[n]*solve(Z2))
                b <- (1-alfa[n])*solve(Z1)%*%m1 + alfa[n]*solve(Z2)%*%m2
                d <- a %*% b
                x <- append(x, d[1], after=length(x))
                y <- append(y, d[2], after=length(y))
            }
        #Create temporary data.frame to make plots easier
        temp_mixture <- data.frame(alfa, (1/2*( dmvnorm(cbind(x,y) , m1, Z1) + dmvnorm(cbind(x,y) , m2, Z2) )))
        names(temp_mixture)[2] <- "pdf"

        #Estimate tolerance ellipses and corresponding beta values for various
        #points along the ridgeline manifold.
        beta_vector_SPA <- c()
        beta_vector_SPB <- c()
        for (j in seq(2, 1000, 1))
            {
                k <- (n_SPA-1)%*%t(c(x[j],y[j])-m1)%*%solve((n_SPA-1)*Z1)%*%(c(x[j],y[j])-m1)
                chi_quantile <- k*(1/(n_SPA-1))*(qchisq(1-gamma, n_SPA-p))
                beta_or_chi_probability_SPA<-pchisq(chi_quantile, p, ncp=p/n_SPA, lower.tail = TRUE, log.p = FALSE)

                k <- (n_SPB-1)%*%t(c(x[j],y[j])-m2)%*%solve((n_SPB-1)*Z2)%*%(c(x[j],y[j])-m2)
                chi_quantile <- k*(1/(n_SPB-1))*(qchisq(1-gamma, n_SPB-p))
                beta_or_chi_probability_SPB<-pchisq(chi_quantile, p, ncp=p/n_SPB, lower.tail = TRUE, log.p = FALSE)

                beta_vector_SPA <- append(beta_vector_SPA, beta_or_chi_probability_SPA, after=length(beta_vector_SPA))
                beta_vector_SPB <- append(beta_vector_SPB, beta_or_chi_probability_SPB, after=length(beta_vector_SPB))
            }

        temp_propA <- data.frame(alfa[2:1000], beta_vector_SPA)
        names(temp_propA)[1] <- "alfa_A"
        names(temp_propA)[2] <- "beta_A"
        temp_propB <- data.frame(alfa[2:1000], beta_vector_SPB)
        names(temp_propB)[1] <- "alfa_B"
        names(temp_propB)[2] <- "beta_B"
        tm1 = t(m1)
        tm2 = t(m2)
        m1 = data.frame(tm1)
        m2 = data.frame(tm2)
        xy = data.frame(x,y)

        #Plots
        mixture_plots[[i]] <- ggplot() +
                                geom_line(aes(alfa, pdf), temp_mixture) +
                                labs( title = paste(tail(sort(subdata$panel), n = 1))) +
                                theme(axis.text = element_text(size = 5),
                                	  axis.title = element_blank(),
                                	  legend.position = "none",
                                	  plot.margin = unit(c(2,1,1,1), "mm"),
                                	  title = element_text( size = 8 ))

        manifold_plots[[i]] <- ggplot(data=subdata, aes(LAMLEN, LAMWID, color = SP, shape = factor(SP))) +
                                geom_point() +
                                geom_point(data = m1, aes(x = X1, y = X2), color = "red", shape = 19, size = 2) +
                                geom_point(data = m2, aes(x = X1, y = X2), color = "red", shape = 19, size = 2) +
                                geom_point(data = xy, aes(x,y), color = "red", shape = 19, size = 2) +
                                labs(x = NULL, y = NULL, title = paste(i)) +
                                theme(axis.text.x = element_text(size = 5),
                                	  axis.text.y = element_text(size = 5),
                                	  legend.position = "none",
                                	  plot.margin = unit(c(2,1,1,1), "mm"),
                                	  title = element_text(size = 10))

        proportion_plots[[i]] <- ggplot() +
                                    ylim(c(0, 1)) +
                                    geom_line(aes(alfa_A, beta_A), temp_propA, color = "#33CCFF" ) +
                                    geom_line(aes(alfa_B, beta_B), temp_propB, color = "000033") +
                                    geom_hline(aes(yintercept=freq_cutoff), linetype = "dashed") +
                                    labs(title=paste(tail(sort(subdata$panel), n=1))) +
                                    theme(axis.text = element_text(size = 5),
                                    	  axis.title = element_blank(),
                                    	  legend.position = "none",
                                    	  plot.margin = unit(c(2,1,1,1), "mm"),
                                    	  title = element_text(size = 8))

        leaf_plots[[i]] <- ggplot(data=subdata, aes(LAMLEN, LAMWID, color = SP, shape = factor(SP))) +
                            geom_point() +
                            labs(x = NULL, y = NULL, title = paste(tail(sort(subdata$panel), n = 1))) +
                            theme(axis.text.x = element_text(size = 5),
                            	  axis.text.y = element_text(size = 5),
                            	  legend.position = "none",
                            	  plot.margin = unit(c(2,1,1,1), "mm"),
                            	  title = element_text(size = 10))

    }

# DISCONTINUITIES IN FLOWER NUMBER

# Are there discontinuities in the number of flowers?

#List to store plot
flowerprop_plots <- list()

#Gamma level of statistical confidence to calculate statistical tolerance regions
q = 0.90

#Frequency cutoff
p = 0.85

for (i in unique(data$species)[unique(data$species) != c("harrisii")])
    {
        subdata = data[data$species %in% c(paste(i), "harrisii"),]

        a <- subset(subdata, species == "harrisii", select = NFLOWERS)
        b <- subset(subdata, species == paste(i), select = NFLOWERS)
        #Some species were not fertile; remove missing data
        a <- na.omit(a)
        b <- na.omit(b)
        mean_a <- mean(a$NFLOWERS)
        mean_b <- mean(b$NFLOWERS)
		sd_a <- sd(a$NFLOWERS)
		sd_b <- sd(b$NFLOWERS)
		n_a <-length(a$NFLOWERS)
		n_b <-length(b$NFLOWERS)

		#Estimate tolerance regions
		A <- 1-( ((qnorm(q))^2)/(2*(n_a-1)) )
		B <- ((qnorm(p))^2)- ( ((qnorm(q))^2)/n_a )
		k1 <- ( qnorm(p) + ((  (qnorm(p)^2)-(A*B) )^0.5) ) / A

		D <- 1-( ((qnorm(q))^2)/(2*(n_b-1)) )
		E <- ((qnorm(p))^2)- ( ((qnorm(q))^2)/n_b )
		k2 <- ( qnorm(p) + ((  (qnorm(p)^2)-(D*E) )^0.5) ) / D

		x_a <- (mean_a + k1*(var(a)^0.5))
		x_b <- (mean_b - k2*(var(b)^0.5))

		flowerprop_plots[[i]] = ggplot(subdata, aes(x = factor(SP), y = NFLOWERS, fill = SP)) +
									geom_boxplot() +
									guides(fill = TRUE) +
									labs(title = paste(tail(sort(subdata$panel), n = 1))) +
									theme(axis.text.y = element_text(size=5),
										  axis.title.y =element_blank(),
										  axis.title.x =element_blank(),
										  axis.text.x = element_blank(),
										  axis.ticks.x = element_blank(),
										  legend.position = "none",
										  plot.margin = unit(c(2,1,1,1), "mm"),
										  title = element_text(size = 10)) +
									geom_hline(yintercept=c(x_a[1], x_b[1]),
											   linetype=c("dashed", "dotted"),
											   color=c("#33CCFF", "#000033"))
	}



-------------------------
#################
#               #
#    FIGURES    #
#               #
#################

=====
## FIGURE 1

args_leaf_plots = c(leaf_plots,
		    list(nrow = 7,
			 ncol = 6,
			 left = "Lamina Width (mm)",
			 bottom = "Lamina Length (mm)"))

do.call(grid.arrange, args_leaf_plots)

=====
## FIGURE 2
args_proportion_plots = c(proportion_plots,
			  list(nrow = 7,
			       ncol = 6,
			       left = "Proportion within tolerance region",
			       bottom = "Ridgeline manifold"))

do.call(grid.arrange, args_proportion_plots)


=====
## FIGURE 3

ggplot(data, aes( x = species, y = NFLOWERS )) +
	   geom_boxplot( na.rm = T ) +
	   theme( axis.text.x  = element_text( angle = 90, vjust = 0.5, size = 10 ),
		 axis.title.y = element_text(size = 15 )) +
	   scale_x_discrete( limits = species_names2,
			    labels = species_names3) +
       labs( x="", y = "Number of flowers" )

=====
## FIGURE 4

args_flowerprop_plots = c(flowerprop_plots,
			  list(nrow = 7,
			       ncol = 6,
			       left = "Number of flowers",
			       bottom = "Species"))

do.call(grid.arrange, args_flowerprop_plots)

=====
## FIGURE 5

ggplot() +
	geom_polygon( data = southam,
		     aes( long, lat, group = group ),
		     fill = NA, color = "grey30" ) +
	geom_point( aes( long, lat),
		   subset( data, !(species %in% out), size = 10, alpha = 0.8, na.rm = T ), 
		   col = "grey65", 
		   size = 1.5 ) +
	geom_point( aes( long, lat, color = species, shape = species ),
		   subset( data, species %in% out, size = 10, na.rm = T ), 
		   size = 3, 
		   alpha = 0.8 )  +
	    scale_color_manual( name = "Species",
                        limits = out,
                        labels = out_names,
                        values = c( "#009E73", "#56B4E9",
                                    "#D55E00", "#000000",
                                    "#E69F00", "#CC79A7",
                                    "#0072B2", "#F0E442" )) +
    scale_shape_manual( name = "Species",
                        limits = out,
                        labels = out_names,
                        values =c( 15, 10,
				  19, 12,
				  17, 18,
				  3, 25)) +
	labs( x = "longitude (degrees)", y="latitude (degrees)" ) +
	theme( panel.background = element_rect( fill = NA ),
		   panel.grid.major = element_line( colour = "grey90" ),
		   axis.title = element_text( size = 10 ),
		   legend.position = "bottom",
		   legend.key = element_rect(fill = NA, color = NA) )


elevp =
	ggplot( data, aes( species, elev, fill = species )) +
    geom_boxplot( aes( species, elev ), 
		 subset( data, species %in% out), 
		 alpha = 0.8, na.rm = T ) +
    theme( axis.text.x  = element_blank(),
    	   axis.ticks.x = element_blank(),
    	   axis.title.y = element_text(size=10)) +
    scale_x_discrete( limits = out, labels = out_names ) +
    scale_fill_manual(values = c( "#56B4E9", "#D55E00",
                                  "#009E73", "#000000",
                                  "#E69F00", "#CC79A7",
                                  "#0072B2", "#F0E442" )) +
    labs(x="", y = "Elevation") +
    #scale_fill_brewer( palette = "Set1" ) +
    theme( legend.position='none',
    	  plot.background = element_rect( colour = "black" ))


vp1 = viewport( width = 0.4,
                height = 0.37,
                x = 0.24,
               	y = 0.4)

print( elevp, vp = vp1 )


=====
## FIGURE 6

# Hand drawing by Barbara Alongi

=====
## FIGURE 7

southam_map =
	ggplot( ) +
	geom_polygon( data = southam,
		     aes(x = long, y = lat, group = group),
		     fill = NA,
		     color = "grey30" ) +
    geom_polygon( data = bolivia,
		 aes(x=long, y = lat, group = group),
		 fill = "gray60") +
    theme( panel.background = element_rect(fill = NA),
		   panel.grid.major = element_line(colour = "grey90"),
		   axis.text = element_blank(),
		   axis.title = element_blank(),
		   axis.ticks = element_blank())

bolivia_map =
    ggplot() +
    geom_polygon( data = bolivia,
		 aes( x = long, y = lat, group = group),
		 col = "gray40",
		 fill = NA) +
	geom_rect( data = zoom_box,
		  aes( xmin = xmin,
		      xmax = xmax,
		      ymin = ymin,
		      ymax = ymax ),
		  fill = "gray70") +
	theme( panel.background = element_rect(fill = NA),
		   panel.grid.major = element_line(colour = "grey90"),
		   axis.text = element_blank(),
		   axis.title = element_blank(),
		   axis.ticks = element_blank())

eh_map2 =
	ggmap(map2) +
    geom_point( aes( x = long, y = lat ),
    			data = subdata,
    			color="black",
    			size = 7 ) +
    geom_point( aes( x = long, y = lat ),
    			data = subdata,
    			color= "white",
    			size = 5.5 ) +
    labs( x = "longitude (degrees)",
    	  y ="latitude (degrees)",
    	  size = 20 ) +
    theme( axis.text = element_text(size = 12),
    	   axis.title = element_text(size=16)) +
    scale_bar( lon = -65.5,
	      lat = -21.45,
	      distance_lon = 50,
	      distance_lat = 5,
	      distance_legend = 10,
	      dist_unit = "km",
	      orientation = FALSE,
	      legend_size = 4,
	      legend_colour = "white" )


vp2 = viewport( width = 0.32,
                height = 0.32,
                x = 0.26,
                y = 0.81)

vp3 = viewport( width = 0.25,
                height = 0.32,
                x = 0.855,
                y = 0.81)
eh_map2
print( southam_map, vp = vp2 )
print( bolivia_map, vp = vp3 )
