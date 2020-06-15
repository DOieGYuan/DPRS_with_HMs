# shell
keggbr="br08001"
echo "level	entry	name" > ${keggbr%}.formatted.txt
cat ${keggbr%}.txt | sed 's/^A/A\t/g' | sed 's/^B */B\t/g' | sed 's/^C */C\t/g' | sed 's/^D */D\t/g' | sed 's/^E */E\t/g' | sed 's/  /\t/g' >> ${keggbr%}.formatted.txt
conda deactivate
R
# R
library(tidyverse)
file = "br08001"
br <- read_tsv(paste(file,".formatted.txt",sep=""))
la <- as.numeric(row.names(br)[br$level=="A"])
la[length(la)+1] <- length(br$entry)+1
br <- mutate(br, L1 = NA)
f <- 1
for(i in la){
  name <- br$entry[f];
  if(i==f){
    f <- i
  }
  else{
    while(f<i){
    br$L1[f] <- name
    f <- f+1
    }
  }
}
lb <- as.numeric(row.names(br)[br$level=="B"])
lb[length(lb)+1] <- length(br$entry)+1
br <- mutate(br, L2 = NA)
f <- 1
for(i in lb){
  name <- br$entry[f];
  if(i==f){
    f <- i
  }
  else{
    while(f<i){
    br$L2[f] <- name
    f <- f+1
    }
  }
}
lc <- as.numeric(row.names(br)[br$level=="C"])
lc[length(lc)+1] <- length(br$entry)+1
br <- mutate(br, L3 = NA)
f <- 1
for(i in lc){
  name <- br$entry[f];
  if(i==f){
    f <- i
  }
  else{
    while(f<i){
    br$L3[f] <- name
    f <- f+1
    }
  }
}
br1 <- filter(br, level == "D" | level == "E") %>% select(-level) %>% mutate(brite=file)
write_tsv(br1, paste(file,".aligned.formatted.txt",sep = ""))
# remove level A,B,C
