traitstats.rbdsummary<-function(Treatment,Replication,DataFile)
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
  GV <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  gv<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      gv=c((TrMSS[i,j]-ErMSS[i,j])/r)
      GV[i,j]<-gv
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
  PV <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  pv<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      pv=c(GV[i,j]+ErMSS[i,j])
      PV[i,j]<-pv
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
  EV <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  ev<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      ev=c(PV[i,j]-GV[i,j])
      EV[i,j]<-ev
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
  GCV <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  gcv<-c()

  for(i in 1:i)
  {
    C<-sqrt(GV[i,i])
    gcv<-c((C/TM[i,i])*100)
    GCV[i,i]<-gcv
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
  GCVCAT <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:i)
  {
    if(GCV[i,i]=="NaN")
    {
      GCVCAT[i,i]<-"ERROR"
    }
    else
    {
      if(GCV[i,i]>20)
      {
        GCVCAT[i,i]<-"HIGH"
      }
      else
      {
        if(10 <  GCV[i,i] && GCV[i,i] < 20)
        {
          GCVCAT[i,i]<-"MODERATE"
        }
        else if(GCV[i,i]<10)
        {
          GCVCAT[i,i]<-"LOW"
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
  PCV <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  pcv<-c()

  for(i in 1:i)
  {
    C<-sqrt(PV[i,i])
    pcv<-c((C/TM[i,i])*100)
    PCV[i,i]<-pcv
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
  PCVCAT <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:i)
  {
    if(PCV[i,i]=="NaN")
    {
      PCVCAT[i,i]<-"ERROR"
    }
    else
    {
      if(PCV[i,i]>20)
      {
        PCVCAT[i,i]<-"HIGH"
      }
      else
      {
        if(10 <  PCV[i,i] && PCV[i,i] < 20)
        {
          PCVCAT[i,i]<-"MODERATE"
        }
        else if(PCV[i,i]<10)
        {
          PCVCAT[i,i]<-"LOW"
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
  h2 <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  h<-c()
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      h=c((GV[i,j]/PV[i,j])*100)
      h2[i,j]<-h
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
  H2CAT <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:i)
  {
    for(j in 1:j)
    {
      if(h2[i,j]=="NaN")
      {
        H2CAT[i,i]<-"ERROR"
      }
      else
      {
        if(h2[i,j] > 60)
        {
          H2CAT[i,j]<-"HIGH"
        }
        else
        {
          if(30 < h2[i,j] && h2[i,j] < 60)
          {
            H2CAT[i,j]<-"MODERATE"
          }
          else if(h2[i,j]<30)
          {
            H2CAT[i,j]<-"LOW"
          }
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
  GA <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  ga<-c()
  for(i in 1:i)
  {
    ga=c((GV[i,i]/PV[i,i])*(2.06)*sqrt(PV[i,i]))
    GA[i,i]<-ga
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
  GACAT <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:i)
  {
    if(GA[i,i]=="NaN")
    {
      GACAT[i,i]<-"ERROR"
    }
    else
    {
      if(GA[i,i] > 20)
      {
        GACAT[i,i]<-"HIGH"
      }
      else
      {
        if(10 <  GA[i,i] && GA[i,i] < 20)
        {
          GACAT[i,i]<-"MODERATE"
        }
        else if(GA[i,i] < 10)
        {
          GACAT[i,i]<-"LOW"
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
  GAM <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  gam<-c()
  for(i in 1:i)
  {
    gam=c((GA[i,i]/TM[i,i])*100)
    GAM[i,i]<-gam
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
  GAMCAT <- matrix( nrow = matrow,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:i)
  {
    if(GAM[i,i]=="NaN")
    {
      GAMCAT[i,i]<-"ERROR"
    }
    else
    {
      if(GAM[i,i] > 20)
      {
        GAMCAT[i,i]<-"HIGH"
      }
      else
      {
        if(10 <  GAM[i,i] && GAM[i,i] < 20)
        {
          GAMCAT[i,i]<-"MODERATE"
        }
        else if(GAM[i,i]<10)
        {
          GAMCAT[i,i]<-"LOW"
        }
      }
    }
  }
  rownames<-c("GCV","PCV","h2","GA","GAM")
  colnames<-c()

  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  GenPar <- matrix( nrow = 5,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:5)
  {
    for(j in 1:j)
    {
      GenPar[1,j]<-paste(round(GCV[j,j],digits = 4),GCVCAT[j,j],sep = ":")
      GenPar[2,j]<-paste(round(PCV[j,j],digits = 4),PCVCAT[j,j],sep = ":")
      GenPar[3,j]<-paste(round(h2[j,j],digits = 4),H2CAT[j,j],sep = ":")
      GenPar[4,j]<-paste(round(GA[j,j],digits = 4),GACAT[j,j],sep = ":")
      GenPar[5,j]<-paste(round(GAM[j,j],digits = 4),GAMCAT[j,j],sep = ":")
    }
  }
  GenPar<-noquote(GenPar)
  rownames<-c("Rep","Trt","MStrt","MSrep","MSerror")
  colnames<-c()

  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  tANOVA <- matrix( nrow = 5,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  j=length((colnames))
  for(i in 1:5)
  {
    for(j in 1:j)
    {

      tANOVA[1,j]<-r
      tANOVA[2,j]<-t
      tANOVA[3,j]<-paste(round(TrMSS[j,j],digits = 4),TSIGN[j,j])
      tANOVA[4,j]<-paste(round(RMSS[j,j],digits = 4),RSIGN[j,j])
      tANOVA[5,j]<-round(ErMSS[j,j],digits = 4)
    }
  }
  tANOVA<-noquote(tANOVA)
  rownames<-c("Trait Mean","Min","Max","SE(m)","CV","CD")
  colnames<-c()
  for(i in 1:matrow)
  {
    colnames <- c(colnames,colnms[i+2])
  }
  tStatsTable <- matrix( nrow = 6,ncol = matrow, byrow = TRUE, dimnames = list(rownames, colnames))
  for(i in 1:6)
  {
    for(j in 1:j)
    {
      tStatsTable[1,j]<-round(TM[j,j],digits=4)
      tStatsTable[2,j]<-round(MIN[j,j],digits=4)
      tStatsTable[3,j]<-round(MAX[j,j],digits=4)
      tStatsTable[4,j]<-round(SEm[j,j],digits=4)
      tStatsTable[5,j]<-round(CV[j,j],digits=4)
      tStatsTable[6,j]<-round((SEd[j,j]*TINV),digits=4)
    }
  }
  tStatsTable

  datarows<-nrow(df_dataset)
  length(colnames)
  dataepr<-c()
  for(j in colrangestart:colcount)
  {
    for(i in 1:datarows)
    {
      dataepr<-as.numeric(append(dataepr,df_dataset[i,j]))
    }
  }
  start = 1
  tastart=1
  tstatstart=1
  tstatend=6
  taend=5
  end = datarows
  plotlist<-list()
  for(i in 1:length(colnames))
  {
    dataplot<-dataepr[start:end]
    plotter<-densityplot(dataplot,col=c("blue"),cex=0.5,lwd=2,xlab="",ylab="",main=colnames[i])
    plotanova <- tableGrob(tANOVA[tastart:taend],rows=rownames(tANOVA))
    plotstats <- tableGrob(tStatsTable[tstatstart:tstatend],rows=rownames(tStatsTable))
    plotgenpar <- tableGrob(GenPar[tastart:taend],rows=rownames(GenPar))

    gp1<-grid.arrange(plotter,ncol=1,nrow=1)
    gp2<-grid.arrange(plotanova,plotstats,ncol=2,nrow=1)
    gp3<-grid.arrange(plotgenpar,ncol=1,nrow=1)
    gp23<-grid.arrange(gp2,gp3,ncol=1,nrow=2)
    g1<-grid.arrange(gp1,gp23,widths=c(2,2),bottom=textGrob("P-value <0.001 *** P-value <0.01 ** P-value <0.05 * ns-non-significant",gp = gpar(fontface = 3, fontsize = 9)))

    plotlist<-list.append(plotlist,g1)
    g1<-list()
    tstatstart=tstatend+1
    tstatend=tstatend+6
    tastart=taend+1
    taend=taend+5
    start = end + 1
    end = end+datarows
  }
  pdfval<-length(plotlist)/2
  looppdfval<-ceiling(pdfval)
  pagestart=1
  pageend=2
  pdfname<-'output.pdf'
  pdf(pdfname)
  for(i in 1:looppdfval)
  {
    if(pageend>length(plotlist))
    {
      pageend=length(plotlist)
    }
    do.call(grid.arrange, c(plotlist[pagestart:pageend],ncol=1,nrow=2))
    pagestart=pagestart+2
    pageend= pageend +2
  }
  dev.off()

  TT1<-grid.text("R:TraitStats Package",gp=gpar(fontsize=40,col="Dark Green"))
  TT2<-grid.text("RCBD Data Analysis Report",gp=gpar(fontsize=30,col="Blue"))
  TT3<-grid.text("Abbreviation:\n ANOVA: Rep-Number of Replication; Trt-Number of Treatment; MStrt-Mean Sum of Squares of Treatment;\n MSrep-Mean Sum of Squares of Replication; MSerror-Mean Sum of Squares of Error.\n DESCRIPTIVE STATISTICS: Trait Mean-Grand Mean of the Trait; Min-Minimum Value; Max-Maximum Value;\n SE(m)-Standard Error of Mean; CV-Coefficient of Variation (%); CD-Critical Difference at 95%.\n GENETIC PARAMETER: GCV-Genotypic Coefficient of variation(%);PCV-Phenotypic Coefficient of variation(%)\n h2-Broad-sense heritability; GA-Genetic Advance; GAM-Genetic Advance percent Mean.",gp=gpar(fontsize=9,col="Black"))
  TT4<-grid.text("Citation:\nNitesh, S.D., Parashuram Patroti and Shilpa Parashuram. (2020).\n TraitStats:Statistical Data Analysis for Randomized Block Design Experiments. R package version 1.0.0",gp=gpar(fontsize=9,col="Black"))

  frontpage<-'fpage.pdf'
  pdf(frontpage)
  do.call(grid.arrange, list(TT1,TT2,TT3,TT4,nrow=4,ncol=1))
  dev.off()

  pdf_combine(c(frontpage, pdfname), output = "TraitStatsRCBD.pdf")
  unlink(frontpage)
  file.remove(pdfname)
  file.show('TraitStatsRCBD.pdf')
}
