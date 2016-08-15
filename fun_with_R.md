### Target file: sam file
```
SRR2426363.195	99	gi|254160123|ref|NC_012967.1|	4017624	60	251M	=	4017900	504	TCATGGAGCTGAAGTCAGTCGAAGATACCAGCTGGCTGCAACTGTTTATTAAAAACACAGCACTGTGCAAACACGAAAGTGGACGTATACGGTGTGACGCCTGCCCGGTGCCGGAAGGTTAATTGATGGGGTCAGCGCAAGCGAAGCTCCTGATCGAAGCCCCGGTAAACGGCGGCCGTAACTATAACGGTCCTAAGGTAGCGAAATTCCTTGTCGGGTAAGTTCCGACCTGCACGAATGGCGTAATGATG	ABABBFFFF?FFCDAFGFEGFEFFGHHDFHGEHHBEBFGGHHHH5FGH5GHFBHG21GHHHHFFHFGHHFHHFGEEG?E@@DFGEEGFHHHG2GEEFF?1EGCCG/EG@EEBECGEDGGFBHDHHG3DHAFD?GGCE@D/@<EGGGFHBGBFFFFGGHFGHGGGGACFFHGGGGGFDFDFDF/BFFFFF;ADCFEFFF9FF/@BB;AFF/BBF:FFFFFA/9/BFFFBBC?FFFFEF?EFBBDA=9EFFFF	XT:A:U	NM:i:3	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:3	XO:i:0	XG:i:0	MD:Z:2G10A135T101
SRR2426363.195	147	gi|254160123|ref|NC_012967.1|	4017900	60	228M	=	4017624	-504	AGTGAAATTGAACTCGCTGTGAAGATGCAGTGTACCCGCGGCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACTGAACATTGAGCCTTGATGTGTAGGATAGGTGGGAGGCTTTGAAGCGTGGACGCCAGTCTGCGTGGAGCCGACCTTGAAATCCCACCCTTTAATGTTTGATGTTCTAACGTGGACCCGTGATCCGGGTTGCGGCCAGTGTCTGGT	;/B/9///FFE=;--E://9/F/9//F9/EFF99-@@G@EGGGFFFBFGGFFE:CA?GEFHG<0=FDD=FGG0HHGBHFDGGD<GF<FHBFAGHHFGFFGBHFFGHCGF1G?G??/?1F0GHGHCFCCCHCB/EE//GFG/?GGGG1/GGCGEFHHHGFG1BEF/CFGGHHGHFGGADF1BFHHHFFFFGGEFGA//E//FHGEEFEAGEAA1EF11BFFFACAAAA?	XT:A:U	NM:i:5	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:5	XO:i:0	XG:i:0	MD:Z:123A4T16A18A51A11
SRR2426363.262	81	gi|254160123|ref|NC_012967.1|	3379269	37	153M	=	3379269	-153	TTGTATTACTCCTCCGACTTACTCAGCGCCGCCGACGAAGTCCAGATTCTGGCCTTCTTTCAGGGTGACGTAAGCTTTTTTCCAGTCGCTACGACGACCGATACGCTGTCCGTGACGTTTAACTTTCCCTTTAACGACCAGGGTGTTAACGAC	GFFGGFGGDGGGGGGHHEHGGHGF?GGGGGGGGEFFFHHHHHHGHHHHHHGFHHHHHHHHGHHF/GHGGFHHHHGGHGFFHGGGGGGGGFGGGFGGGGGHDDEGHGGGGGGGCHHGHGHHHHHHFHGHHHFGGGGGFGGGGCFFFCCBCCCBB	XT:A:U	NM:i:2	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:2	XO:i:0	XG:i:0	MD:Z:33A101T17
```

```
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
#---  summarise(): reduce each group to a smaller number of summary statistics
#---  arrange(): re-order the rows
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
    ## return(x)
  }) %>%
  `rownames<-`(unlist(QNAME))
```  


<hr />
Let's break it down!
<hr />

```{R}

##==== Step 1: split each line by tab and return a tabular data structure ====
myData %>% tail(-1) %>% ## remove the first row because it contains field names
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
##==== Step 2: set column names ====
  (function(x){
    colnames(x) = c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                    'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL')
    return(x)
  }) %>%
```

```{R}
##==== Step:3 convert strings to numbers ====
  mutate(FLAG = as.numeric(FLAG),
         pos = as.numeric(POS),
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
