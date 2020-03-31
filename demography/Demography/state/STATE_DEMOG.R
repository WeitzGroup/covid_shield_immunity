
# dAta is from
# https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-detail.html


FILE = "PEP_2018_PEPSYASEX_with_ann.csv"
header<- scan(FILE,nlines = 1, what = character())
head2 = strsplit(header,",")
AA<- read.csv(FILE, stringsAsFactors = FALSE, header=FALSE, skip=2, )
colnames(AA) = head2[[1]]


#data is in 1 year age groups up to 85 (which is 85+)



#est4 ...  are summaries
#est7 ... are age based descriptions
#sex0 ... everyone  (this is what we want!)
#sex1 ... denoted male
#sex2 ... denoted female



#Want columns of est72018sex0_ageX  where X is age (0,...,85)
StringStart = "est72018sex0_age"
Ages = seq(0,85,by=1)
StringsOfInterest = paste(StringStart,Ages,sep="")
StringsOfInterest[86] = paste(StringsOfInterest[86],"plus",sep="")


rNAMES = AA[,3] #states 

NEWDATA = matrix(0,length(rNAMES),length(Ages))

#subset data
for(ii in 1:length(Ages)){
IND = which(head2[[1]] == StringsOfInterest[ii])
NEWDATA[,ii] = AA[,IND]
}


#recast into decade (10 year age groups)
#0-9,10-19,20-29,30-39,40-49,50-59,60-69,70-79,80+
MinAge = c(0,10,20,30,40,50,60,70,80)
NEWDATA2 = matrix(0,length(rNAMES),length(MinAge))


for(aa in 1:(length(MinAge)-1)){
NEWDATA2[,aa] =  rowSums(NEWDATA[,( (aa-1)*10+1 ):( (aa-1)*10+10 )])
}

NEWDATA2[,9] = rowSums(NEWDATA[,81:86])

POPTOTALS = rowSums(NEWDATA2)


POPFRACS = 0.*NEWDATA2
for(ii in 1:length(rNAMES)){
POPFRACS[ii,] = NEWDATA2[ii,] /POPTOTALS[ii]
}


STATE_DEMOG = cbind(POPFRACS,POPTOTAL)

rownames(STATE_DEMOG) = rNAMES
colnames(STATE_DEMOG) = c("0to9","10to19","20to29","30to39","40to49","50to59","60to69","70to79","80+","POPTOTAL")

write.csv(STATE_DEMOG,"STATE_DEMOG.csv")





