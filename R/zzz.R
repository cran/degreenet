.onAttach <- function(lib, pkg){
  desc <- packageDescription("degreenet")
  copylist <- "Copyright (c) 2013, Mark S. Handcock, University of California -- Los Angeles"
  sm <- paste("\n",desc$Package,": version ", desc$Version, ', created on ', desc$Date, '\n',copylist,"\n",
       'Based on "statnet" project software (statnet.org).\n',
       'For license and citation information see statnet.org/attribution\n',
       'or type citation("',desc$Package,'").\n', sep="")
  if(!is.null(sm)) packageStartupMessage(sm)
}
