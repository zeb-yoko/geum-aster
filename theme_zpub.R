##Publication Theme adapted from Jill's
##black and white##
theme_zpub <- function (base_size = 12, base_family = "") 
{
	theme_classic(base_size = base_size, base_family = base_family) %+replace% 
		theme(axis.text = element_text(size = rel(0.8)), 
				axis.ticks = element_line(colour = "black"), 
				legend.key = element_rect(colour = "gray43"), 
				panel.border = element_rect(fill = NA, colour = "black"))
}

##set colors	#Prairie   #GL ALvar  #MB Alvar##
##col.esa <- c("#18563E", "#82BE42", "#FFC423")
##col.grscl <- c("gray80", "black", "gray43")