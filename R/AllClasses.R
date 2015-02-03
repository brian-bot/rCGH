# CLASS DEFINITION

setClass('cghObj', representation(	info = 'character',			# a vector of char containing.... specify what!
													cnSet = 'data.frame',							# a data.frame containing the CN probes values
													param = 'vector',								  # a vector of mixted values (boolean/numeric) containing the parameters used by CBS algorithm.
													segTable = 'data.frame',					# a data.frame containing the segment values. Availabe after the segmenation step.
                          byGene = 'data.frame',					  # a data.frame containing the genes values. Availabe after the segmenation step.
                          probesDensity = 'ANY',						# a xyplot (of class 'trellis') representing the probe density as a gaussian mixture.
													gProfile = 'ANY'),									# a xyplot (of class 'trellis') representing the genomic profile.
         prototype=prototype(info = c(),
                             cnSet = data.frame(),
                             param = NA,
                             segTable = data.frame(),
                             byGene = data.frame(),
                             probesDensity = NULL,
                             gProfile = NULL
#                              probesDensity = xyplot(c(0, 1)~c(0, 1), type  ='n'),
#                              gProfile = xyplot(c(0, 1)~c(0, 1), type  ='n')
                             )
						)

# SHOW METHOD FOR THIS CLASS
setMethod('show',
          signature = 'cghObj',
          definition = function(object){
            d <- dim(object@cnSet)
            infoTab = data.frame(info = getInfo(object))
            cat('\nInstance of class', class(object), '\n\n')
            cat('CNSet with', d[1], 'probes and', d[2], 'columns\n\n')
            cat('Array information:\n\n')
            for(i in 1:nrow(infoTab)){
              item = rownames(infoTab)[i]
              value = as.character(infoTab[i,1])
              cat('\t', item, ":", value, '\n')
            }
            cat('\n')
            cat('Use getInfo(object) to get array information\n')
            cat('Use getCNset(object) to get the CGH matrix\n')
            cat('Use getParam(object) to get segmentation parameters\n')
            cat('Use getSegTable(object) to get segmentation table\n')
            cat('Use getProfile(object) to display the genomic profile\n\n')
          }
)
