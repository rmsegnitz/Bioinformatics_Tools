
#### Deep save ggplot objects to overwrite without conflict ####

  gg_deepCopy<-function(x){
    
    if(all(class(x) %in% c("gg", "ggplot"))){
      return(unserialize(serialize(x, NULL)))
    }else{
      print("Error: object must be gg or ggplot class object")
    }
    
  }