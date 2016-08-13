```{R}
library(plyr)
library(dplyr)
library(stringr)
rm(list=ls())
myData = readLines('https://raw.githubusercontent.com/MingChen0919/KBS_workshop_Michigan_2016/master/SRR2426363.mapped.sam')
myData = readLines('SRR2426363.mapped.sam')
```

```{R}
## Exploring
#---  select(): focus on a subset of variables
#---  filter(): focus on a subset of rows
#---  mutate(): add new columns
#---  summarise(): reduce each group to a smaller number of summary statistics
#---  arrange(): re-order the rows


## CIGAR
#------  Op BAM Description
#------  M 0 alignment match (can be a sequence match or mismatch)
#------  I 1 insertion to the reference
#------  D 2 deletion from the reference
#------  N 3 skipped region from the reference
#------  S 4 soft clipping (clipped sequences present in SEQ)
#------  H 5 hard clipping (clipped sequences NOT present in SEQ)
#------  P 6 padding (silent deletion from padded reference)
#------  = 7 sequence match
#------  X 8 sequence mismatch
```

```{R}
## read file by lines
myData %>% 
  (function(x){
    ## wrap the str_split_fixed() function so it returns a character
    str_split_fixed_2 = function(string, pattern, n){
      x = str_split_fixed(string, pattern, n)
      return(as.character(x))
    }
    ldply(x, str_split_fixed_2, '\t', 20)[, 1:11]
  }) %>% 
  (function(x){
    colnames(x) = c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                    'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL')
    return(x)
  }) %>% 
  mutate(FLAG = as.numeric(FLAG),
         pos = as.numeric(POS),
         MAPQ   = as.numeric(MAPQ),
         PNEXT  = as.numeric(PNEXT),
         TLEN   = as.numeric(TLEN)) %>%
  (function(x){
    filter(x, MAPQ > 30 & MAPQ < 55) %>%
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>% nrow() 
    return(x)
  }) %>% 
  (function(x){
    arrange(x, MAPQ) %>% 
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>% head() 
    return(x)
  }) %>%
  (function(x){
    filter(x, MAPQ > 30) %>% 
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>%  
      arrange(POS, MAPQ) %>% head() 
    return(x)
  }) %>% 
  (function(x){
    select(x, CIGAR) %>% 
      laply(grep, pattern='[I]') %>% 
      slice(.data = x) %>%
      (function(x){
        x ->> df_with_Insertion
        select(x, CIGAR) ->> CIGAR_with_Insertion
        select(x, QNAME) ->> QNAME 
        return(select(x, CIGAR))
      })
  }) %>% 
  llply(str_split, 'I') %>%
  `[[`('CIGAR') %>%
  llply(head, -1) %>%
  llply(gsub, pattern='[A-Z]', replacement='+') %>%
  (function(x){
    eval_string = function(var1) eval(parse(text=var1))  ## define a function
    eval_string_plus = function(var2) laply(var2, eval_string) ## define another function
    llply(x, eval_string_plus)
  })  %>% 
  llply(cumsum) %>%
  (function(x){
    maxN = max(lengths(x))
    newX = vector('numeric', length = maxN)
    llply(x, `[`, 1:maxN)
  }) %>% as.data.frame() %>% t() %>%
  (function(x){
    `colnames<-`(x, paste0('insert_', 1:ncol(x)))
  }) %>%
  `rownames<-`(QNAME)
```  




```{R}
## read file by lines
myData %>% 
  (function(x){
    ## wrap the str_split_fixed() function so it returns a character
    str_split_fixed_2 = function(string, pattern, n){
      x = str_split_fixed(string, pattern, n)
      return(as.character(x))
    }
    ldply(x, str_split_fixed_2, '\t', 20)[, 1:11]
  }) %>%
```

```{R}
  (function(x){
    colnames(x) = c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                    'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL')
    return(x)
  }) %>%
```

```{R}
  mutate(FLAG = as.numeric(FLAG),
         pos = as.numeric(POS),
         MAPQ   = as.numeric(MAPQ),
         PNEXT  = as.numeric(PNEXT),
         TLEN   = as.numeric(TLEN)) %>%
```

```{R}
  (function(x){
    filter(x, MAPQ > 30 & MAPQ < 55) %>%
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>% nrow() 
    return(x)
  }) %>%
```

```{R}
  (function(x){
    arrange(x, MAPQ) %>% 
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>% head() 
    return(x)
  }) %>%
```  

```{R}  
  (function(x){
    filter(x, MAPQ > 30) %>% 
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>%  
      arrange(POS, MAPQ) %>% head() 
    return(x)
  }) %>%
```  

```{R}  
  (function(x){
    select(x, CIGAR) %>% 
      laply(grep, pattern='[I]') %>% 
      slice(.data = x) %>%
      (function(x){
        x ->> df_with_Insertion
        select(x, CIGAR) ->> CIGAR_with_Insertion
        select(x, QNAME) ->> QNAME 
        return(select(x, CIGAR))
      })
  }) %>%
```

```{R}
  llply(str_split, 'I') %>%
```

```{R}
  `[[`('CIGAR') %>%
```

```{R}
  llply(head, -1) %>%
```

```{R}
  llply(gsub, pattern='[A-Z]', replacement='+') %>%
```

```{R}
  (function(x){
    eval_string = function(var1) eval(parse(text=var1))  ## define a function
    eval_string_plus = function(var2) laply(var2, eval_string) ## define another function
    llply(x, eval_string_plus)
  })  %>% 
```

```{R}
  llply(cumsum) %>%
```
```{R}
  (function(x){
    maxN = max(lengths(x))
    newX = vector('numeric', length = maxN)
    llply(x, `[`, 1:maxN)
  }) %>% as.data.frame() %>% t() %>%
```

```{R}
  (function(x){
    `colnames<-`(x, paste0('insert_', 1:ncol(x)))
  }) %>%
  `rownames<-`(QNAME)
```  
