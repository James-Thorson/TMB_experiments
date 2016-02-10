
MakeADFun_windows_debug = function( data, parameters, random, cpp_name, dir=getwd(), ... ){
  # Get and save inputs
  Other_inputs = list(...)
  All_inputs = list( "data"=data, "parameters"=parameters, "random"=random, "Other_inputs"=Other_inputs )
  save( All_inputs, file="All_inputs.RData")

  # Write file to source
  sink( paste0(cpp_name,".R") )
  cat("
    library( TMB )
    dyn.load( dynlib('",cpp_name,"') )
    setwd('",dir,"')
    load( 'All_inputs.RData' )
    Obj = MakeADFun(data=All_inputs[['data']], parameters=All_inputs[['parameters']], random=All_inputs[['random']], All_inputs[['Other_inputs']])
    save( Obj, file='Obj.RData')
  ",fill=FALSE,sep="")
  sink()

  #
  gdbsource(paste0(cpp_name,".R"))
  #if( !is(Try, "try-error")
  #save( x, file=paste0(x,".RData") )
  #load( file='Obj.RData')
  #return( Obj )
}
