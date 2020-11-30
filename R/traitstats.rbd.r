traitstats.rbd<-function(Treatment,Replication,DataFile)
{
  df_dataset<-DataFile
  treatment<-as.factor(Treatment)
  replication<-as.factor(Replication)

  colnms<-colnames(df_dataset)
  colcount=length(colnms)
  colrangestart=3

  t=length(levels(treatment))
  r=length(levels(replication))

  tlevels=levels(treatment)
  rlevels=levels(replication)

  TFT95<-qf(0.95,(t-1),(r-1)*(t-1))
  TFT99<-qf(0.99,(t-1),(r-1)*(t-1))
  TFT99.99<-qf(0.999,(t-1),(r-1)*(t-1))

  TFR95<-qf(0.95,(r-1),(r-1)*(t-1))
  TFR99<-qf(0.99,(r-1),(r-1)*(t-1))
  TFR99.99<-qf(0.999,(r-1),(r-1)*(t-1))

  TINV<-abs(qt(0.05/2,(r-1)*(t-1)))

  sumvector<-c()
  for (val in 3:colcount)
  {
    sumval<-sum(df_dataset[val])
    sumvector<-c(sumvector,sumval)
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RTS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  tempvar=colcount-2
  temprc=length(sumvector)
  for(i in 1:tempvar)
  {
    varstart=i
    for(j in varstart:tempvar)
    {
      if(i==j)
      {
        valTS=0
        {
          valTS= valTS+(sumvector[i]*sumvector[i])
        }
        RTS[i,j]<-valTS
      }
      else
      {
        valTS=0
        valTS= valTS+(sumvector[i]*sumvector[j])
        RTS[i,j]<-valTS
        RTS[j,i]<-valTS
      }

    }
  }
  totvals=t*r
  rtss<-c()
  i=colrangestart
  while(i<=colcount)
  {
    fsum<-c()
    for(k in 1:totvals)
    {
      if(i<=colcount)
      {
        val1=df_dataset[k,i]*df_dataset[k,i]
        fsum <-as.numeric(c(fsum,val1))
      }
    }
    tss=sum(as.numeric(fsum))
    rtss <-c(rtss,tss)
    tempi=i
    while(tempi<=colcount)
    {
      fsum<-c()
      for(k in 1:totvals)
      {
        if(tempi<=colcount-1)
        {
          val1=df_dataset[k,i]*df_dataset[k,tempi+1]
          fsum <-as.numeric(c(fsum,val1))
        }
      }
      tempi=tempi+1
      TSS=sum(fsum)

      rtss <-c(rtss,TSS)
    }
    i=i+1
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RTSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  tempvar=colcount-2
  countchar=1
  for(i in 1:tempvar)
  {
    varstart=i
    for(j in varstart:tempvar)
    {
      if(rtss[countchar]==0)
        countchar=countchar+1

      {
        val=rtss[countchar]
        RTSS[i,j]<-val
        RTSS[j,i]<-val
        countchar=countchar+1
      }
    }
  }
  repsum<-c()
  rcount=length(rlevels)
  for(col in 3:colcount)
  {
    repcolsum<-c()
    for (valRS in 1:rcount)
    {
      sumval<-sum(df_dataset[which(df_dataset[,2]==rlevels[valRS]),col])
      repcolsum<-c(repcolsum,sumval)
    }
    repsum<-c(repsum,repcolsum)
  }
  rownames<-c()
  colnames<-c()
  repcount=length(rlevels)
  for(i in 1:repcount)
  {
    rownames <- c(rownames,rlevels[i])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RRSS1 <- matrix( nrow = r,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))

  incrementor=1
  cc=colcount-2
  for(i in 1:cc)
  {
    for(j in 1:r)
    {
      RRSS1[j,i]=repsum[incrementor]
      incrementor=incrementor+1
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RRSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  tempvar=colcount-2
  for(i in 1:tempvar)
  {
    varstart=i
    for(j in varstart:tempvar)
    {
      if(i==j)
      {
        valRS=0
        for(k in 1:r)
        {
          valRS=valRS+(RRSS1[k,j]*RRSS1[k,j])
        }
        RRSS[i,j]<-valRS
      }
      else
      {
        valRS=0
        for(k in 1:r)
        {
          valRS=valRS+(RRSS1[k,i]*RRSS1[k,j])
        }
        RRSS[i,j]<-valRS
        RRSS[j,i]<-valRS
      }
    }
  }
  trsum<-c()
  trtcount=length(tlevels)
  for(col in 3:colcount)
  {
    trcolsum<-c()
    for (val in 1:trtcount)
    {
      sumval<-sum(df_dataset[which(df_dataset[,1]==tlevels[val]),col])
      trcolsum<-c(trcolsum,sumval)
    }
    trsum<-c(trsum,trcolsum)
  }
  rownames<-c()
  colnames<-c()
  for(i in 1:trtcount)
  {
    rownames <- c(rownames,tlevels[i])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RTSS1 <- matrix( nrow = t,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  incrementor=1
  cc=colcount-2
  for(i in 1:cc)
  {
    for(j in 1:t)
    {
      RTSS1[j,i]=trsum[incrementor]
      incrementor=incrementor+1
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RTrSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))

  tempvar=colcount-2
  temprc=length(tlevels)

  for(i in 1:tempvar)
  {
    varstart=i
    for(j in varstart:tempvar)
    {
      if(i==j)
      {
        valTS=0
        for(k in 1:temprc)
        {
          valTS= valTS+(RTSS1[k,i]*RTSS1[k,i])

        }
        RTrSS[i,j]<-valTS
      }
      else
      {
        valTS=0
        for(k in 1:temprc)
        {
          valTS= valTS+(RTSS1[k,j]*RTSS1[k,i])

        }
        RTrSS[i,j]<-valTS
        RTrSS[j,i]<-valTS
      }

    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  CFM <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  CF<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      CF=c(RTS[i,j]/(t*r))
      CFM[i,j]<-CF
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  TSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  TS<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      TS=c(RTSS[i,j]-CFM[i,j])
      TSS[i,j]<-TS
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  TrSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  TrS<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      TrS=c((RTrSS[i,j]/(r))-CFM[i,j])
      TrSS[i,j]<-TrS
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  RS<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      RS=c((RRSS[i,j]/(t))-CFM[i,j])
      RSS[i,j]<-RS
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  ErSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  ErS<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      ErS=c(TSS[i,j]-TrSS[i,j]-RSS[i,j])
      ErSS[i,j]<-ErS
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  TrMSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  TrMS<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      TrMS=c(TrSS[i,j]/(t-1))
      TrMSS[i,j]<-TrMS
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RMSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  RMS<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      RMS=c(RSS[i,j]/(r-1))
      RMSS[i,j]<-RMS
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  ErMSS <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  ErMS<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      A<-(t-1)*(r-1)
      ErMS=c(ErSS[i,j]/A)
      ErMSS[i,j]<-ErMS
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  CalFT <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  calft<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      calft=c(TrMSS[i,j]/ErMSS[i,j])
      CalFT[i,j]<-calft
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  CalFR <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  calfr<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      calfr=c(RMSS[i,j]/ErMSS[i,j])
      CalFR[i,j]<-calfr
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  TSIGN <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      if(TFT99.99 < CalFT[i,j])
      {
        TSIGN[i,j]<-"***"
      }
      else
      {
        if(TFT99<CalFT[i,j])
        {
          TSIGN[i,j]<-"**"
        }
        else if(TFT95<CalFT[i,j])
        {
          TSIGN[i,j]<-"*"
        }
        else
        {
          TSIGN[i,j]<-"ns"
        }
      }
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  RSIGN <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      if(TFR99.99 < CalFR[i,j])
      {
        RSIGN[i,j]<-"***"
      }
      else
      {
        if(TFR99<CalFR[i,j])
        {
          RSIGN[i,j]<-"**"
        }
        else if(TFR95<CalFR[i,j])
        {
          RSIGN[i,j]<-"*"
        }
        else
        {
          RSIGN[i,j]<-"ns"
        }
      }
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  SEm <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  sem<-c()
  for(i in 1:i)
  {
    E=(ErMSS[i,i]/r)
    sem=c(sqrt(E))
    SEm[i,i]<-sem
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  SEd <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  sed<-c()
  for(i in 1:i)
  {
    E=((2*ErMSS[i,i])/r)
    sed=c(sqrt(E))
    SEd[i,i]<-sed
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  TM <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  traittotal=sumvector
  tm<-c()
  for(i in 1:i)
  {
    for(j in i:i)
    {
      tm<-traittotal[i]/(t*r)
      TM[i,j]<-c(tm)
    }
  }
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      tm<-c(TM[i,i]+TM[j,j])/2
      TM[i,j]<-tm
      TM[j,i]<-tm
    }
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  MIN <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  mind<-c()
  for(i in 1:i)
  {
    mind<-c(min(df_dataset[i+2]))
    MIN[i,i]<-mind
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  MAX <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  maxd<-c()
  for(i in 1:i)
  {
    maxd<-c(max(df_dataset[i+2]))
    MAX[i,i]<-maxd
  }
  matrow=colcount-2
  rownames<-c()
  colnames<-c()
  for(i in 1:matrow)
  {
    rownames <- c(rownames,colnms[i+2])
  }
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  CV <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  cv<-c()
  for(i in 1:i)
  {
      B<-sqrt(ErMSS[i,i])
      cv<-c((B/TM[i,i])*100)
      CV[i,i]<-cv
  }
  rownames<-c("Number of Replication","Number of Treatment","Treatment Mean Sum of Square","Replication Mean Sum of Square","Error Mean Sum of 		Square","Trait Mean","Minimum","Maximum","CV","StdError(m)","StdError(d)")
  colnames<-c()

  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  ANOVA <- matrix( nrow = 11,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
j=length((colnames))
           for(i in 1:11)
           {
             for(j in 1:j)
             {
               ANOVA[1,j]<-r
               ANOVA[2,j]<-t
               ANOVA[3,j]<-paste(round(TrMSS[j,j],digits = 4),TSIGN[j,j])
               ANOVA[4,j]<-paste(round(RMSS[j,j],digits = 4),RSIGN[j,j])
               ANOVA[5,j]<-round(ErMSS[j,j],digits=4)
               ANOVA[6,j]<-round(TM[j,j],digits=4)
               ANOVA[7,j]<-round(MIN[j,j],digits=4)
               ANOVA[8,j]<-round(MAX[j,j],digits=4)
               ANOVA[9,j]<-round(CV[j,j],digits=4)
               ANOVA[10,j]<-round(SEm[j,j],digits=4)
               ANOVA[11,j]<-round(SEd[j,j],digits=4)
             }
           }
           ANOVA<-noquote(ANOVA)
           ANOVA
           return(ANOVA)
}
