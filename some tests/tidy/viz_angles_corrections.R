library(magrittr)
library(purrr)

viz_angles_corrections<-function(ref, ...){
  list(...) -> corrections
  corrections %>% map(~class(.x[[1]])) %>% equals("list") -> parametric
  ggplot(data.frame(
    k = corrections[parametric] %>% map(seq_along) %>% unlist,#à changer
    correction = corrections[parametric] %>% imap(~.y %>% rep(length(.x))) %>% unlist,#à changer
    angle = corrections[parametric] %>% map(~.x %>% map(~.x %>% map(~.x[[1]]))) %>% unlist,
    batch = corrections[parametric] %>% map(~.x %>% map(names)) %>% unlist
  ))
}