sig_gg_network_manual_pies<- # overwrite with deep copy & reorder layers
  unserialize(serialize(move_layers(sig_gg_network_manual_pies, "GeomTextRepel", position = "top"), NULL))

#### Deep save ggplot objects to overwrite without conflict ####

  gg_deepCopy<-function(x){
    
    if(all(class(x) %in% c("gg", "ggplot"))){
      return(unserialize(serialize(x)))
    }else{
      print("Error: object must be gg or ggplot class object")
    }
    
  }