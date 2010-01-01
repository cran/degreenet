#  File R/zzz.R
#  Part of the "statnet" package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California-Los Angeles
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Martina Morris, University of Washington
# Copyright 2007 statnet Development Team
#
# We have invested a lot of time and effort in creating 'statnet',
# for use by other researchers. We require that the attributions
# are retained with each function.
######################################################################
######################################################################
#
# .First.lib is run when the package is loaded.
#

.First.lib <- function(lib, pkg){
    library.dynam("degreenet", pkg, lib)
    DESCpath <- file.path(system.file(package="degreenet"), "DESCRIPTION")
    info <- read.dcf(DESCpath)
    cat('\ndegreenet:', info[,"Title"], 
        '\nVersion', info[,"Version"], 'created on', info[,"Date"], '\n') 
    cat(paste("copyright (c) 2003, Mark S. Handcock, University of California-Los Angeles\n",
"                    David R. Hunter, Penn State University\n",
"                    Carter T. Butts, University of California-Irvine\n",
"                    Steven M. Goodreau, University of Washington\n",
"                    Pavel N. Krivitsky, Carnegie Mellon University\n",
"                    Martina Morris, University of Washington\n",sep=""))
    cat('Type help(package="degreenet") to get started.\n\n')
    cat('Based on "statnet" project software (http://statnet.org).\n',
        'For license and citation information see http://statnet.org/attribution\n',
        'or type citation("degreenet").\n')
}

.Last.lib <- function(libpath){
  library.dynam.unload("degreenet",libpath)
}
