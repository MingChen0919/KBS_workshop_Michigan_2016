### Target file: sam file
```
SRR2426363.615	163	gi|254160123|ref|NC_012967.1|	3356120	10	71M1D8M1I125M23S	=	3356375	504	CTTGCACCCTCCGTATTACCGCGGCTGCTGGCACGGAGTTAGCCGGTGCTTCTTCTGCGGGTAACGTCAATTGCTGCGGTTATTAGCCACAACACCTTCCTCCCCGCTGAAAGTACTTTACAACCCGAAGGCCTTCTTCATACACGCGGCATGGCTGCATCAGGCTTGCGCCCATTGTGCAATATTCCCCACTGCTGCCTCCCGTCGTAGTCAGTACTGTGTCTGAGT	AAABBFF43AAAEGFGGGGGFGGFCGGFHHHGHHGGGG?GHHHHC?FC?G5GHGGGGHGF/EEHH?GGGCFHHHHGFGGGG?GBFGEF3?C?<GE?BGGHHHH/E?CCD/DGFDGHHGHBGHGFD.<A@C?.C<<0DGFFFHFFGEGF?D-AFGGGEEGGFF9FB?BBEA?D?.//9/9BFBFFFFFF/;BAFEFBBF.;9.9.A.;-99AB/;/9;/:99BBB//:9	XT:A:M	NM:i:13	SM:i:10	AM:i:10	XM:i:11	XO:i:2	XG:i:2	MD:Z:71^G0A2A0A0A7A1T0T0T1C0T0C111
SRR2426363.698	89	gi|254160123|ref|NC_012967.1|	4015323	0	168M1I33M	=	4015323	0	TGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTACCTTAAAGAAGCGTACTTTGCAGTGCTCACACAGATTGTCTGATGAAAAACGAGCAGTAAAACCTCTACAGGCTTGTAG	EA;0C00GFC:00C:/DGGCFD0D0DDC?GGCC.<<--/@@??1?11?1GCGDCDHFDFB@@1B<?@FB12DCFAEF>F11<>1DHE>>AEHFFGF1EF0EEE/>/B1B1HG1E>EFBFFFGFGFHHDGBEFAA/2GFGAEGDF1GG1GFB2ABGF1HGFCB1F2GGF0F0GEEAEDFG1DF1B1FABF3FFFAAC1?AAAA	XT:A:R	NM:i:2	SM:i:0	AM:i:0	X0:i:3	X1:i:0	XM:i:1	XO:i:1	XG:i:1	MD:Z:172T28	XA:Z:gi|254160123|ref|NC_012967.1|,-228028,168M1I33M,2;gi|254160123|ref|NC_012967.1|,+3355048,28M1I173M,2;
SRR2426363.1034	163	gi|254160123|ref|NC_012967.1|	3912080	29	39M1I137M	=	3912459	618	CTCTTTAGGTATTCCTTCGAACAAGATGCAAGAATAGACAAAAATGACAGCCCTTCTACGAGTGATTAGCCTGGTCGTGATTAGCGTGGTGGTGATTATTATCCCACCGTGCGGGGCTGCACTTGGACGAGGAAAGGCTTAGAAATCAAGCCTTAACGAACTAAGACCCCCGCACCG	ABBBBFFBFFFFGGGGGGGCGGGHGHFFFFCHHGHHHHFHFFAFGHFFFGEAEHHHHHGFDEEHHHHFHHHGFGGHFEGHHHHGF0B????EEFHHGGGHHFGHHHHGGHHEE?EC?DGFBGFHHEFFCGG//CGC0CGEB11=1>GG1>GGBGFB0DC.-.F0;CFHEGCCGCA@9	XT:A:U	NM:i:4	SM:i:29	AM:i:29	X0:i:1	X1:i:0	XM:i:3	XO:i:1	XG:i:1	MD:Z:9C24A107G33
```


![CIGAR](https://raw.githubusercontent.com/MingChen0919/KBS_workshop_Michigan_2016/master/Screen%20Shot%202016-08-16%20at%204.08.57%20PM.png)


## Goals
### Positions of all insertion/deletion
```{R}
insert_1	insert_2	insert_3	insert_4
SRR2426363.615	81	NA	NA	NA
SRR2426363.698	169	NA	NA	NA
SRR2426363.1034	40	NA	NA	NA
SRR2426363.4217	164	176	NA	NA
```
### Data manipulation skills with __plyr/dplyr__
```{R}
#---  **ply(): a*ply, l*ply, d*ply
#---  select(): focus on a subset of variables
#---  filter(): focus on a subset of rows
#---  mutate(): add new columns
#---  arrange(): re-order the rows
#---  [[: subsetting
#---  ->: right assignment operator
#---  ->>: assign through parents environments
#---  %>%
#---  workflow
```


<hr />
__Let's do it!__
<hr />

## Load packages
```{R}
p = c('plyr', 'dplyr', 'stringr')
install.packages(p)

library(plyr)
library(dplyr)
library(stringr)
rm(list=ls())
```

## Get data
### Option 1: read data remotely
```{R}
myData = readLines('https://raw.githubusercontent.com/MingChen0919/KBS_workshop_Michigan_2016/master/SRR2426363.mapped.sam')
```

### Option 2: download data and then read data locally
```{R}
curl -O https://raw.githubusercontent.com/MingChen0919/KBS_workshop_Michigan_2016/master/SRR2426363.mapped.sam
myData = readLines('SRR2426363.mapped.sam')
```

### Do everything in one single pipe line.

* Pipeline starts with your data
* Input -> Output: No intermediate variables generated
* Pipeline is extensible

```{R}
myData %>% tail(-1) %>% ## remove the first row because it contains field names
  (function(x){
    colnames(x) = c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                    'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL')
    return(x)
  }) %>% 
  (function(x){
    ## wrap the str_split_fixed() function so it returns a character
    str_split_fixed_2 = function(string, pattern, n){
      x = str_split_fixed(string, pattern, n)
      return(as.character(x))
    }
    ldply(x, str_split_fixed_2, '\t', 20)[, 1:11]
  }) %>% 
  mutate(FLAG = as.numeric(FLAG),
         POS = as.numeric(POS),
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
    ## return(x)
  }) %>%
  `rownames<-`(unlist(QNAME))
```  


<hr />
Let's break it down!
<hr />

```{R}

##==== Step 1:  ====
myData %>% tail(-1) %>% ## remove the first row because it contains field names
  (function(x){ ## set column names ====
    colnames(x) = c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                    'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL')
    return(x)
  }) %>%
 
```

```{R}
==== Step 2: split each line by tab and return a tabular data structure ===
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
##==== Step:3 convert strings to numbers ====
  mutate(FLAG = as.numeric(FLAG),
         POS = as.numeric(POS),
         MAPQ   = as.numeric(MAPQ),
         PNEXT  = as.numeric(PNEXT),
         TLEN   = as.numeric(TLEN)) %>%
```

```{R}
##==== Exploring  ====
  (function(x){
    filter(x, MAPQ > 30 & MAPQ < 55) %>%
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>% nrow() 
    return(x)
  }) %>%
```

```{R}
##==== Exploring ====
  (function(x){
    arrange(x, MAPQ) %>% 
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>% head() 
    return(x)
  }) %>%
```  

```{R}
##==== Exploring ====
  (function(x){
    filter(x, MAPQ > 30) %>% 
      select(QNAME, RNAME, POS, MAPQ, CIGAR) %>%  
      arrange(POS, MAPQ) %>% head() 
    return(x)
  }) %>%
```  

```{R}
##==== Step 4: select rows that have insertions ====
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
##==== Step 5: split CIGAR by insertion letter 'I'  ====
  llply(str_split, 'I') %>%
```

```{R}
##==== Step 6: Extract data from a nested list ====
  `[[`('CIGAR') %>%
```

```{R}
##==== Step 7: delete the last element ====
  llply(head, -1) %>%
```

```{R}
##==== Step 8: build calculation expression by replace remaining CIGAR letters with '+' ====
  llply(gsub, pattern='[A-Z]', replacement='+') %>%
```

```{R}
##==== Step 9: evalate the expression ====
  (function(x){
    eval_string = function(var1) eval(parse(text=var1))  ## define a function
    eval_string_plus = function(var2) laply(var2, eval_string) ## define another function
    llply(x, eval_string_plus)
  })  %>% 
```

```{R}
##==== Step 10: get positions ====
  llply(cumsum) %>%
  (function(x){
    maxN = max(lengths(x))
    ## newX = vector('numeric', length = maxN)
    llply(x, `[`, 1:maxN)
  }) %>% 
```  

```{R} 
##==== Step 10: Convert list to data frame ====
  as.data.frame() %>% t() %>%
```

```{R}
##==== Step 11: column names and row names ====
  (function(x){
    `colnames<-`(x, paste0('insert_', 1:ncol(x)))
  }) %>%
  `rownames<-`(unlist(QNAME))
```  
