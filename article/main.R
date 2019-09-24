library(magrittr)
library(purrr)
library(eigenangles)
library(ggplot2)

#Mouse
experiments <- load_experiments('Mouse Expression Atlas/batch data/')
experiments %<>% remove_isolated_experiments('organism_part')
experiments %<>% merge_experiments
experiments %>% save(file='Mouse Expression Atlas/integrated data/uncorrected.Rdata')

experiments %>% correct_batch_effect(~organism_part, 'ComBat') %>% 
  save(file='Mouse Expression Atlas/integrated data/ComBat_corrected.Rdata')
1:20 %>% map(~experiments %>% correct_batch_effect(~organism_part, 'RUV', k=.x)) %>% 
  save(file='Mouse Expression Atlas/integrated data/RUV_corrected_k=1:20.Rdata')
1:20 %>% map(~experiments %>% correct_batch_effect(~organism_part, 'MNN', k=.x)) %>% 
  save(file='Mouse Expression Atlas/integrated data/MNN_corrected_k=1:20.Rdata')
experiments %>% correct_batch_effect(~organism_part, 'BMC') %>% 
  save(file='Mouse Expression Atlas/integrated data/BMC_corrected.Rdata')

apply_eigenangles(
  uncorrected = get(load('Mouse Expression Atlas/integrated data/uncorrected.Rdata')),
  ComBat = get(load('Mouse Expression Atlas/integrated data/ComBat_corrected.Rdata')),
  RUVs = get(load('Mouse Expression Atlas/integrated data/RUVs_corrected_k=1:20.Rdata')),
  MNN = get(load('Mouse Expression Atlas/integrated data/MNN_corrected_k=1:20.Rdata')),
  BMC = get(load('Mouse Expression Atlas/integrated data/BMC_corrected.Rdata')),
  group = 'organism_part'
) -> ea

ea %>% save(file='Mouse Expression Atlas/integrated data/eigenangles.Rdata')

ea %>% plot+
  scale_colour_manual(values=c('orange','red','green','blue','black'))
ea %>% tanmean %>% plot+
  scale_colour_manual(values=c('orange','red','green','blue','black'))

#Human
experiments <- load_experiments('Human Expression Atlas/batch data/')
experiments %<>% remove_isolated_experiments('organism_part')
experiments %<>% merge_experiments
experiments %>% save(file='Human Expression Atlas/integrated data/uncorrected.Rdata')

experiments<-get(load('Human Expression Atlas/integrated data/uncorrected.Rdata'))
experiments %>% correct_batch_effect(~organism_part, 'ComBat') %>% 
  save(file='Human Expression Atlas/integrated data/ComBat_corrected.Rdata')
1:5 %>% map(~experiments %>% correct_batch_effect(~organism_part, 'RUV', k=.x)) %>% 
  save(file='Human Expression Atlas/integrated data/RUVs_corrected_k=1:5.Rdata')
1:5 %>% map(~experiments %>% correct_batch_effect(~organism_part, 'MNN', k=.x)) %>% 
  save(file='Human Expression Atlas/integrated data/MNN_corrected_k=1:5.Rdata')
experiments %>% correct_batch_effect(~organism_part, 'BMC') %>% 
  save(file='Human Expression Atlas/integrated data/BMC_corrected.Rdata')

apply_eigenangles(
  uncorrected = get(load('Human Expression Atlas/integrated data/uncorrected.Rdata')),
  ComBat = get(load('Human Expression Atlas/integrated data/ComBat_corrected.Rdata')),
  RUVs = get(load('Human Expression Atlas/integrated data/RUVs_corrected_k=1:5.Rdata')),
  MNN = get(load('Human Expression Atlas/integrated data/MNN_corrected_k=1:5.Rdata')),
  BMC = get(load('Human Expression Atlas/integrated data/BMC_corrected.Rdata')),
  group = 'organism_part'
) -> ea

ea %>% save(file='Human Expression Atlas/integrated data/eigenangles.Rdata')

ea %>% plot+
  scale_colour_manual(values=c('orange','red','green','blue','black'))
ea %>% tanmean %>% plot+
  scale_colour_manual(values=c('orange','red','green','blue','black'))



#gPCA
#Mouse
get(load('Mouse Expression Atlas/integrated data/uncorrected.Rdata')) %>% do_gPCA->gun
get(load('Mouse Expression Atlas/integrated data/ComBat_corrected.Rdata')) %>% do_gPCA->gcb
get(load('Mouse Expression Atlas/integrated data/RUVs_corrected_k=1:20.Rdata')) %>% map(do_gPCA)->gruv
get(load('Mouse Expression Atlas/integrated data/MNN_corrected_k=1:20.Rdata')) %>% map(do_gPCA)->gmnn
get(load('Mouse Expression Atlas/integrated data/BMC_corrected.Rdata')) %>% do_gPCA->gbmc

gun %>% save(file='Mouse Expression Atlas/gPCA/uncorrected.Rdata')
gcb %>% save(file='Mouse Expression Atlas/gPCA/ComBat.Rdata')
gruv %>% save(file='Mouse Expression Atlas/gPCA/RUVs_k=1:20.Rdata')
gmnn %>% save(file='Mouse Expression Atlas/gPCA/MNN_k=1:20.Rdata')
gbmc %>% save(file='Mouse Expression Atlas/gPCA/BMC.Rdata')

gun <- get(load('Mouse Expression Atlas/gPCA/uncorrected.Rdata'))
gcb <- get(load('Mouse Expression Atlas/gPCA/ComBat.Rdata'))
gruv <- get(load('Mouse Expression Atlas/gPCA/RUVs_k=1:20.Rdata'))
gmnn <- get(load('Mouse Expression Atlas/gPCA/MNN_k=1:20.Rdata'))
gbmc <- get(load('Mouse Expression Atlas/gPCA/BMC.Rdata'))

gun %>% replace(c('upca','gpca','batch'),NULL) %>% save(file='Mouse Expression Atlas/gPCA_statistics/uncorrected.Rdata')
gcb %>% replace(c('upca','gpca','batch'),NULL) %>% save(file='Mouse Expression Atlas/gPCA_statistics/ComBat.Rdata')
gruv %>% map(~.x %>% replace(c('upca','gpca','batch'),NULL)) %>% save(file='Mouse Expression Atlas/gPCA_statistics/RUVs_k=1:20.Rdata')
gmnn %>% map(~.x %>% replace(c('upca','gpca','batch'),NULL)) %>% save(file='Mouse Expression Atlas/gPCA_statistics/MNN_k=1:20.Rdata')
gbmc %>% replace(c('upca','gpca','batch'),NULL) %>% save(file='Mouse Expression Atlas/gPCA_statistics/BMC.Rdata')

gun <- get(load('Mouse Expression Atlas/gPCA_statistics/uncorrected.Rdata'))
gcb <- get(load('Mouse Expression Atlas/gPCA_statistics/ComBat.Rdata'))
gruv <- get(load('Mouse Expression Atlas/gPCA_statistics/RUVs_k=1:20.Rdata'))
gmnn <- get(load('Mouse Expression Atlas/gPCA_statistics/MNN_k=1:20.Rdata'))
gbmc <- get(load('Mouse Expression Atlas/gPCA_statistics/BMC.Rdata'))

B=7
K=20
ggplot(
  data.frame(
    algorithm=c(
      'BMC',
      'ComBat',
      'MNN' %>% rep(K),
      'RUV' %>% rep(K),
      'uncorrected'
    ),
    k=c(
      NA,
      NA,
      1:K,
      1:K,
      NA
    ),
    delta=c(
      gbmc$delta[1],
      gcb$delta[1],
      gmnn %>% map(~.x$delta[1]) %>% unlist,
      gruv %>% map(~.x$delta[1]) %>% unlist,
      gun$delta[1]
    )
  )
)+aes(x=k,y=delta,colour=algorithm)+geom_point()


ggplot(
  data.frame(
    algorithm=c(
      'BMC' %>% rep(B),
      'MNN' %>% rep(B*K),
      'RUV' %>% rep(B*K),
      'ComBat' %>% rep(B),
      'uncorrected' %>% rep(B)
    ),
    k=c(
      NA %>% rep(B),
      1:K %>% rep(each=B),
      1:K %>% rep(each=B),
      NA %>% rep(B),
      NA %>% rep(B)
    ),
    delta=c(
      gbmc$delta,
      gmnn %>% map(~.x$delta) %>% unlist,
      gruv %>% map(~.x$delta) %>% unlist,
      gcb$delta,
      gun$delta
    ),
    dim=c(
      1:B,
      1:B %>% rep(K),
      1:B %>% rep(K),
      1:B,
      1:B
    )
  )
)+aes(x=dim,y=delta,colour=algorithm,alpha=k)+geom_point()+
  scale_colour_manual(values=c('orange','red','green','blue','black'))


#Human
get(load('Human Expression Atlas/integrated data/uncorrected.Rdata')) %>% do_gPCA->gun
get(load('Human Expression Atlas/integrated data/ComBat_corrected.Rdata')) %>% do_gPCA->gcb
get(load('Human Expression Atlas/integrated data/RUVs_corrected_k=1:5.Rdata')) %>% map(do_gPCA)->gruv
get(load('Human Expression Atlas/integrated data/MNN_corrected_k=1:5.Rdata')) %>% map(do_gPCA)->gmnn
get(load('Human Expression Atlas/integrated data/BMC_corrected.Rdata')) %>% do_gPCA->gbmc

gun %>% save(file='Human Expression Atlas/gPCA/uncorrected.Rdata')
gcb %>% save(file='Human Expression Atlas/gPCA/ComBat.Rdata')
gruv %>% save(file='Human Expression Atlas/gPCA/RUVs_k=1:5.Rdata')
gmnn %>% save(file='Human Expression Atlas/gPCA/MNN_k=1:5.Rdata')
gbmc %>% save(file='Human Expression Atlas/gPCA/BMC.Rdata')

gun <- get(load('Human Expression Atlas/gPCA/uncorrected.Rdata'))
gcb <- get(load('Human Expression Atlas/gPCA/ComBat.Rdata'))
gruv <- get(load('Human Expression Atlas/gPCA/RUVs_k=1:5.Rdata'))
gmnn <- get(load('Human Expression Atlas/gPCA/MNN_k=1:5.Rdata'))
gbmc <- get(load('Human Expression Atlas/gPCA/BMC.Rdata'))

gun %>% replace(c('upca','gpca','batch'),NULL) %>% save(file='Human Expression Atlas/gPCA_statistics/uncorrected.Rdata')
gcb %>% replace(c('upca','gpca','batch'),NULL) %>% save(file='Human Expression Atlas/gPCA_statistics/ComBat.Rdata')
gruv %>% map(~.x %>% replace(c('upca','gpca','batch'),NULL)) %>% save(file='Human Expression Atlas/gPCA_statistics/RUVs_k=1:5.Rdata')
gmnn %>% map(~.x %>% replace(c('upca','gpca','batch'),NULL)) %>% save(file='Human Expression Atlas/gPCA_statistics/MNN_k=1:5.Rdata')
gbmc %>% replace(c('upca','gpca','batch'),NULL) %>% save(file='Human Expression Atlas/gPCA_statistics/BMC.Rdata')

gun <- get(load('Human Expression Atlas/gPCA_statistics/uncorrected.Rdata'))
gcb <- get(load('Human Expression Atlas/gPCA_statistics/ComBat.Rdata'))
gruv <- get(load('Human Expression Atlas/gPCA_statistics/RUVs_k=1:5.Rdata'))
gmnn <- get(load('Human Expression Atlas/gPCA_statistics/MNN_k=1:5.Rdata'))
gbmc <- get(load('Human Expression Atlas/gPCA_statistics/BMC.Rdata'))

B=5
K=5
ggplot(
  data.frame(
    algorithm=c(
      'BMC',
      'ComBat',
      'MNN' %>% rep(K),
      'RUV' %>% rep(K),
      'uncorrected'
    ),
    k=c(
      NA,
      NA,
      1:K,
      1:K,
      NA
    ),
    delta=c(
      gbmc$delta[1],
      gcb$delta[1],
      gmnn %>% map(~.x$delta[1]) %>% unlist,
      gruv %>% map(~.x$delta[1]) %>% unlist,
      gun$delta[1]
    )
  )
)+aes(x=k,y=delta,colour=algorithm)+geom_point()


ggplot(
  data.frame(
    algorithm=c(
      'BMC' %>% rep(B),
      'MNN' %>% rep(B*K),
      'RUV' %>% rep(B*K),
      'ComBat' %>% rep(B),
      'uncorrected' %>% rep(B)
    ),
    k=c(
      NA %>% rep(B),
      1:K %>% rep(each=B),
      1:K %>% rep(each=B),
      NA %>% rep(B),
      NA %>% rep(B)
    ),
    delta=c(
      gbmc$delta,
      gmnn %>% map(~.x$delta) %>% unlist,
      gruv %>% map(~.x$delta) %>% unlist,
      gcb$delta,
      gun$delta
    ),
    dim=c(
      1:B,
      1:B %>% rep(K),
      1:B %>% rep(K),
      1:B,
      1:B
    )
  )
)+aes(x=dim,y=delta,colour=algorithm,alpha=k)+geom_point()+
  scale_colour_manual(values=c('orange','red','green','blue','black'))
