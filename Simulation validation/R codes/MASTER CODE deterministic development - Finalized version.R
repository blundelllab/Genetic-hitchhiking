
# + Each mutation has an ID (ID) and a growth advantage (s) as a tuple
# ``[(ID1, s1)]``
# + Each clone contains a certain number of mutations, so could be defined as a list of the mutation, fitness effect tuples
# ``[(ID1, s1),(ID2, s2)....]``
# + Each clone also needs an associated size (n) (e.g. 50 cells) -> ``[[(ID1, s1),(ID2, s2)....], n]``
# + Total population is a list of the different clones (with a different set of mutations within it): ``[[[(ID1, s1),(ID2, s2)....], n], [[(ID1, s1),(ID2, s2)....], n]]]``
#this turns off mean fitness update

#import libraries
library(rlist)
library(distr)
library(ggplot2)
library(scales)
library(gridExtra)
library(beepr)


start_time <- Sys.time()

benefit=0.1

#Define function 'Divide'
Divide<- function(clone_info, dt, B0, D0,mean_fitness){
  
  cloneID=clone_info[[1]] # [[ID1, s1],[ID2,s2]...] is the ID of the clone
  clone_size=clone_info[[2]][1] # n
  #     print('clone_size =',clone_size)
  
  Total_fitness_effects=0;
  for (i in 1:length(cloneID) )
  {Total_fitness_effects=Total_fitness_effects+cloneID[[i]][[2]]}
  #fitness of the clone is the sum of fitness conferred by all mutations it contains
  Total_fitness_effects=Total_fitness_effects - mean_fitness                                     
  
  B = B0+Total_fitness_effects #B = birth rate, which is 1 + the overall fitness effect of the mutations within a particular clone
  if (clone_size!=0){
    if(B>0){
      number_of_births <- rpois(1, clone_size*B*dt)#pulls a random number from the Poisson distribution with/
      #mean of clone_size x birth rate x interval of time
      
    }else{number_of_births=0}            #to make sure the mean is positive        
    
    number_of_deaths <- rpois(1, clone_size*D0*dt);
    #     print('number of deaths = ', number_of_deaths)
  }
  else{number_of_births=0;
  number_of_deaths=0; }          #makes sure an extinct clone stays extinct
  
  
  new_clone_size = clone_size + number_of_births - number_of_deaths;
  if(new_clone_size<0){new_clone_size=0}
  #print(paste("the new clone size is ",new_clone_size))
  return (new_clone_size)
}

#Define function 'Divide_dev'
Divide_dev<- function(clone_info, dt, B0, D0, mean_fitness, r){
  
  cloneID=clone_info[[1]] # [[ID1, s1],[ID2,s2]...] is the ID of the clone
  clone_size=clone_info[[2]][1] # n
  
  new_clone_size = clone_size*exp(r*dt);
  
  if(clone_size==0){new_clone_size=0}
  #print(paste("the new clone size is ", new_clone_size))
  return (new_clone_size)
}

#Define function 'mutation_fitness' 
mutation_fitness <- function(ratio_neutral_beneficial,benefit) {
  
  random_number<- runif(1, 0, 1)
  threshold= 1/(1+ratio_neutral_beneficial)
  if(random_number<threshold){
    
    s= benefit

    
  }else{s=0}
  
  return(s)                                  
}


#Define function 'mutate'
mutate<-function(clone_info, dt, u, last_mutation_ID, last_development_mutation_ID,ratio_neutral_beneficial,ratio_unknown_known){
  
  mutations=clone_info[[1]] #[[ID1, s1],[ID2,s2]]
  clone_size=clone_info[[2]]
  
  list_of_new_clones<-list()
  list_of_ben_mut_ID<-list()
  list_of_neu_mut_ID<-list()
  list_of_known_ben_mut_ID<-list()
  
  threshold=1/(ratio_unknown_known+1)
  
  
  number_of_mutations <-rpois(1, clone_size*u*dt) ;
  
  
  if(number_of_mutations!=0){
    
    for (i in 1:number_of_mutations){
      last_mutation_ID = last_mutation_ID+1
      
      new_fitness_effect = mutation_fitness(ratio_neutral_beneficial,benefit)
      
      if(new_fitness_effect==0){
        
        list_of_neu_mut_ID <-list.append(list_of_neu_mut_ID, last_mutation_ID)
        
      }
      else{
        list_of_ben_mut_ID<-list.append(list_of_ben_mut_ID, last_mutation_ID)
        if(runif(1, 0, 1)<threshold){
          list_of_known_ben_mut_ID<-list.append(list_of_known_ben_mut_ID, last_mutation_ID)}
      }
      
      new_mutation<-list(last_mutation_ID, new_fitness_effect)
      new_clone_info <-list(list.append(mutations,new_mutation),1)    #[[[ID1,s1],[ID2,s2],....],n]
      
      list_of_new_clones<-list.append(list_of_new_clones,new_clone_info)
      
    }
  }   
  
  result<-list(list_of_new_clones, last_mutation_ID, number_of_mutations,
               list_of_ben_mut_ID, list_of_neu_mut_ID, list_of_known_ben_mut_ID)     
  return (result)                                     
}


#define function 'update_mean_fitness'
update_mean_fitness<- function(Bio_system)
{        
  
  number_of_clones= length(Bio_system)
  Total_fitness_effects=0;
  total_pop=0;
  weighted_fitness=0;
  for (i in 1:number_of_clones){                #i is the tag of the clone
    number_of_mutations=length(Bio_system[[i]][[1]])  #number of mutations in that clone
    Total_fitness_effects=0
    for (j in 1:number_of_mutations) {
      Total_fitness_effects= Total_fitness_effects+Bio_system[[i]][[1]][[j]][[2]]
    }                                       #sum up the fitness conferred by each mutation in a clone
    
    weighted_fitness=weighted_fitness+Bio_system[[i]][[2]]*Total_fitness_effects
    total_pop=total_pop+Bio_system[[i]][[2]];                 
    
  }
  mean_fitness=weighted_fitness/total_pop;
  
  return (list(mean_fitness, number_of_clones,total_pop))       
}


#Main function


Number_of_Participants=1000
Population_size=10^5                                 #is the starting population size of clone
hitchhiker_cutoff=0.5;
known_cutoff=0.5
ratio_unknown_known=0
dt=0.1;#in units of generation
number_of_generations=70;#number of generations

r=1.2
T_dev=log(Population_size)/(dt*r) # average developmental time

num_non_extinct_lineage=0
lifespan=number_of_generations/dt; #number of runs
Total_mutation_frequency_beneficial<-list()
Total_mutation_frequency_beneficial_first_mutant<-list()
Total_mutation_frequency_beneficial_double_mutant<-list()
Total_mutation_frequency_neutral<-list()
Total_mutation_frequency_neutral_not_hitchhiker<-list()
Total_mutation_frequency_neutral_hitchhiker<-list()
Total_mutation_frequency_neutral_beneficial_first<-list()
Total_mutation_frequency_neutral_beneficial_later<-list()
Total_mutation_frequency_dev<- list()
Total_mutation_frequency_neutral_dev<-list()
Total_mutation_frequency_beneficial_dev<-list()
Total_mutation_frequency_neutral_hitchhiker_developmental_union<-list()
Total_mutation_frequency_neutral_hitchhiker_with_known_ben<-list()
Clonal_histories_per_person<-list()
Entire_experiment<-list()

for (N in 1:Number_of_Participants) {
  
  Current_Population_Size=1  # starts from 1 single cell
  Bio_system <-list(list(list(list(1,0)), Current_Population_Size));
  t=0;
  last_mutation_ID=1;
  last_development_mutation_ID = 0
  B0=1;
  D0=1;
  # Bolton_dev_per_doubling_synonymous_mu_LeeSix = 0.00019816
  u_dev=0.0002*3*r/log(2)   # per time
  u_dev=0.0002*3*r/log(2)/10 
  # Bolton_synonymous_mu 0.000614
  u_ben=0.0006*2*0.05  # 4% nonsyn is above fitness 10%
  u_ben=0.00001
  u_ben_KD=u_ben   # adulthood and development

  u_neu=0.0006*1  # adulthood 
  u= u_ben + u_neu
  ratio_neutral_beneficial=u_neu/u_ben   
  ratio_neutral_beneficial_dev= u_neu/u_ben_KD 
  mean_fitness=0;
  End_Population_Size=0;
  End_of_development=0
  status='development'
  Clonal_histories<-list() 
  beneficial_mut_ID <-list()
  neutral_mut_ID <-list()
  known_beneficial_mut_ID<-list()
  
  Stored_mean_fitness<-list()
  Stored_number_of_clones<-list()
  Bio_system_record<-list()
  mean_fitness_record<-list()
  Number_of_nonextinct_clones_record<-list()
  
  
  for (t in 1:lifespan) {
    
    # calculate current population size
    Current_Population_Size=0
    
    
    for (n in 1:length(Bio_system)) {
      
      Size<- Bio_system[[n]][[2]]
      if(Size!=0){Current_Population_Size=Current_Population_Size+Size}
      
    }
    
    
    

    # makes sure the code knows whether it is still in development
    if (status=='no longer relevant'){
      
      status='no longer relevant'
      
    }else{
      
      if (Current_Population_Size<Population_size){
        
        status='development'
        
      }
      if (Current_Population_Size>Population_size){
        if (status=='development'){
          End_of_development=t*dt
          last_development_mutation_ID = last_mutation_ID
          status='no longer relevant'
          
        }
      }
      
    }
    
    
    
    
   
    
    
    #mutate Bio_system 
    result=list();
    list_of_new_clones=list()
    number_of_mutations=list()
    
    for (m in 1:length(Bio_system)){
      
      if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
        
        if (status=='no longer relevant'){
          result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID, last_development_mutation_ID, ratio_neutral_beneficial,ratio_unknown_known) #consider the first clone
        }else{result[[m]]<- mutate(Bio_system[[m]], dt, u_dev, last_mutation_ID, last_development_mutation_ID, ratio_neutral_beneficial_dev,ratio_unknown_known )}
        
        list_of_new_clones[[m]]<- result[[m]][[1]]
        last_mutation_ID <- result[[m]][[2]]
        number_of_mutations[[m]]<- result[[m]][[3]]
        list_of_ben_mut_ID<-result[[m]][[4]]
        list_of_neu_mut_ID<-result[[m]][[5]]
        list_of_known_ben_mut_ID<-result[[m]][[6]]
        #only call the function 'mutate' once for each clone
        
      }
      if (length(list_of_known_ben_mut_ID)!=0){
        for (i in 1:length(list_of_known_ben_mut_ID)){
          known_beneficial_mut_ID  <-list.append(known_beneficial_mut_ID,list_of_known_ben_mut_ID[[i]])
        }
      }
      
      if (length(list_of_ben_mut_ID)!=0){
        for (i in 1:length(list_of_ben_mut_ID)){
          beneficial_mut_ID  <-list.append(beneficial_mut_ID,list_of_ben_mut_ID[[i]])
        }
      }
      
      if (length(list_of_neu_mut_ID)!=0){
        for (i in 1:length(list_of_neu_mut_ID)){
          neutral_mut_ID <-list.append(neutral_mut_ID, list_of_neu_mut_ID[[i]])
        }
      }
    }
    
    
    
    
    for (m in 1:length(Bio_system)){
      if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones
        
        if (number_of_mutations[[m]]>0){
          for (i in 1:number_of_mutations[[m]]) {
            
            if(length(result[[m]][[1]])!=0){
              
              Bio_system <- list.append(Bio_system, list_of_new_clones[[m]][[i]]) #[[[[ID1, s1],[ID2,s2]],n],[[[ID1,s1],[ID3,s3]],n],[n]]
            }
            
            
          } 
        }
        
      }
    }
    #print("Finished mutating")
    #divide: propagate growth of each clone
    
    
    Clonal_sizes<-list()
    
    for (i in 1:length(Bio_system)) {
      
      if(Bio_system[[i]][[2]]!=0){                           #only divide non-extinct clones
        
        
        if (status=='no longer relevant'){
          Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
        }else{
          Bio_system[[i]][[2]]=Divide_dev(Bio_system[[i]], dt, B0, D0,mean_fitness, r)[[1]]   #update clone size
        }   
        
      }                                                       #consider relative fitness
      Clonal_sizes[[i]]=Bio_system[[i]][[2]];
      
    }
    #print("Finished dividing")
    
    #update mean fitness
    result2<-update_mean_fitness(Bio_system)
    mean_fitness<-result2[[1]]
    Stored_mean_fitness<-list.append(Stored_mean_fitness,mean_fitness)
    
    total_pop<-result2[[2]]
    
    # print(Current_Population_Size>Population_size)
    # print(status=='adulthood')
    
    
    
    
    #print(paste("this is generation ",t))
    
  }
  
  beneficial_mut_ID<-unique(beneficial_mut_ID)
  neutral_mut_ID<-unique(neutral_mut_ID)
  beneficial_dev_mut_ID <- beneficial_mut_ID[beneficial_mut_ID<last_development_mutation_ID]
  neutral_dev_mut_ID <- neutral_mut_ID[neutral_mut_ID<last_development_mutation_ID]
  known_beneficial_mut_ID<-unique(known_beneficial_mut_ID)
  print(paste('Last developmental mutation ID is:', last_development_mutation_ID))
  print(paste('Last mutation ID is:', last_mutation_ID))
  
  #find mutation frequency for this patient
  mutation_frequency_beneficial<-list()                              #all beneficial mutations
  mutation_frequency_beneficial_first_mutant<-list()
  mutation_frequency_beneficial_double_mutant<-list()
  mutation_frequency_neutral<-list()                                #all neutral mutations
  mutation_frequency_neutral_not_hitchhiker<-list()
  mutation_frequency_neutral_hitchhiker<-list()
  mutation_frequency_neutral_beneficial_first<-list()
  mutation_frequency_neutral_beneficial_later<-list()
  mutation_frequency_dev <-list()
  mutation_frequency_beneficial_dev<-list()
  mutation_frequency_neutral_dev<-list()
  mutation_frequency_neutral_hitchhiker_developmental_union<-list()
  mutation_frequency_neutral_hitchhiker_with_known_ben<-list()
  
  # counts end population of the person
  for (n in 1:length(Bio_system)) {
    
    Size<- Bio_system[[n]][[2]]
    if(Size!=0){End_Population_Size=End_Population_Size+Size}
    
  }
  
  if(End_Population_Size!=0){
    num_non_extinct_lineage = num_non_extinct_lineage + 1
  }
  
  # inititalize lists
  for (m in 1:last_mutation_ID){mutation_frequency_beneficial[[m]]=0;
  mutation_frequency_beneficial_first_mutant[[m]]=0;
  mutation_frequency_beneficial_double_mutant[[m]]=0;
  mutation_frequency_neutral[[m]]=0;
  mutation_frequency_neutral_not_hitchhiker[[m]]=0;
  mutation_frequency_neutral_hitchhiker[[m]]=0;
  mutation_frequency_neutral_beneficial_first[[m]]=0;
  mutation_frequency_neutral_beneficial_later[[m]]=0;
  mutation_frequency_dev[[m]]=0
  mutation_frequency_beneficial_dev[[m]]=0
  mutation_frequency_neutral_dev[[m]]=0
  mutation_frequency_neutral_hitchhiker_developmental_union[[m]]=0
  mutation_frequency_neutral_hitchhiker_with_known_ben[[m]]=0
  }    #initialization
  
  for (m in 2:last_mutation_ID){if(last_mutation_ID!=1) {
    
    Size_of_non_hitchhiker_portion=0;
    Size_of_beneficial_first_portion=0;
    Size_of_beneficial_later_portion=0;
    Size_of_with_known_beneficial_mutation=0
    
    for (n in 1:length(Bio_system)) {
      
      CloneID<- Bio_system[[n]][[1]]
      Size<- Bio_system[[n]][[2]]
      
      
      Type_of_neutral_mutationID=0;                           #this can be "Not Hitchhiker", "Beneficial first" or "Beneficial later"
      where_first_beneficial_mutation_comes_in=0;           #to be the position of the 1st beneficial mutation in the n-th clone
      number_of_beneficial_mutations=0;
      clone_has_known_ben_mut=0
      total_fitness=0;
      
      #calculate total fitness of clone
      for (k in 1:length(CloneID)) {
        mutation_ID<- CloneID[[k]][[1]]
        fitness<- CloneID[[k]][[2]]
        
        total_fitness=total_fitness+CloneID[[k]][[2]];
        
        
      }
      
      #determine if the clone is neutral
      if(total_fitness==0){                                #relies on the fact that WT mutation is neutral
        Type_of_neutral_mutationID="Not Hitchhiker"
        #print(paste("not hitchhiker"))
        
      }
      else{
        #if clone not neutral, then find when the first beneficial mutation comes in
        for (k in 1:length(CloneID)) {
          if(CloneID[[k]][[2]]!=0){
            where_first_beneficial_mutation_comes_in=k;         #position of the 1st beneficial mutation
            #print(paste("where_first_beneficial_mutation_comes_in is",where_first_beneficial_mutation_comes_in))
            break
          }
          
        }
        
        
        
        # # #check if clone contains beneficial mutation
        # # this part is EXTREMELY COMPUTATIONAL EXPENSIVE!!!!
        # 
        # for (k in 1:length(CloneID)) {
        # 
        #   if(CloneID[[k]][[1]] %in% known_beneficial_mut_ID){
        #     # clone contains known beneficial mutation
        #     clone_has_known_ben_mut='True'
        #     break
        #   }
        # }
        clone_has_known_ben_mut='True'
        
        #count the number of beneficial mutations in a clone
          for (k in 1:length(CloneID)) {
          if(CloneID[[k]][[2]]!=0){
            number_of_beneficial_mutations=number_of_beneficial_mutations+1;        # number_of_beneficial_mutations in the n-th clone


          }
          
        }
        
      }
      
      #make the lists for mutation frequencies
      
      
      for (k in 1:length(CloneID)) {
        
        if(CloneID[[k]][[1]]==m){
          
          if(m<last_development_mutation_ID){
            # print('developmental!')
            mutation_frequency_dev[[m]]= mutation_frequency_dev[[m]]+Size}
          
          
            if(m %in% beneficial_dev_mut_ID){
            # print(paste('developmental beneficial! size :', Size))
            mutation_frequency_beneficial_dev[[m]]= mutation_frequency_beneficial_dev[[m]]+Size
            }
          
          
            if(m %in% neutral_dev_mut_ID){
            # print('developmental neutral!')
            mutation_frequency_neutral_dev[[m]]= mutation_frequency_neutral_dev[[m]]+Size
            }
          
          
          if(CloneID[[k]][[2]]!=0){#if the mutation is not neutral
            
            
            mutation_frequency_beneficial[[m]]=mutation_frequency_beneficial[[m]]+Size
            
            if(number_of_beneficial_mutations==1)
            {mutation_frequency_beneficial_first_mutant[[m]]=mutation_frequency_beneficial_first_mutant[[m]]+Size}
            if(number_of_beneficial_mutations==2)
            {mutation_frequency_beneficial_double_mutant[[m]]=mutation_frequency_beneficial_double_mutant[[m]]+Size}
            #print(paste("double mutant enters, not necessarily establishes"))
            #this counts the wildtype mutation as well if it is beneficial
            
            #print(paste("mutation ", m,"is beneficial"))
            
          }
          else{
            
            mutation_frequency_neutral[[m]]=mutation_frequency_neutral[[m]]+Size;
            
            
            if(k<where_first_beneficial_mutation_comes_in & where_first_beneficial_mutation_comes_in!=0){ Type_of_neutral_mutationID="Beneficial later"}
            if(k>where_first_beneficial_mutation_comes_in & where_first_beneficial_mutation_comes_in!=0){ Type_of_neutral_mutationID="Beneficial first"}
            if(where_first_beneficial_mutation_comes_in==0){ Type_of_neutral_mutationID="Not Hitchhiker"}
            if(clone_has_known_ben_mut=='True'){Type_of_neutral_mutationID="with known beneficial mutation"}
            
            #print(paste("The neutral mutation ", m," belongs to the type", Type_of_neutral_mutationID, "in clone ",n))
            
            
            
            if( Type_of_neutral_mutationID=="Not Hitchhiker"){
              Size_of_non_hitchhiker_portion=Size_of_non_hitchhiker_portion+Size;}
            if( Type_of_neutral_mutationID=="Beneficial first"){
              Size_of_beneficial_first_portion=Size_of_beneficial_first_portion+Size;}
            if( Type_of_neutral_mutationID=="Beneficial later"){
              Size_of_beneficial_later_portion=Size_of_beneficial_later_portion+Size}
            if(Type_of_neutral_mutationID=="with known beneficial mutation"){
              Size_of_with_known_beneficial_mutation=Size_of_with_known_beneficial_mutation+Size}
          }
          
          break                  
        }        #no need to look at the rest of the CloneID since a mutation only occurs once
        
        
        
      }# end for hitchhiker type
      
    }# end of going through each clone
    # print('end of going through each clone')
    
    if(mutation_frequency_neutral[[m]]!=0){  #it is possible that that a mutation is found in two clones, a neutral clone and a hitchhiker clone
      ratio1 <-Size_of_non_hitchhiker_portion/mutation_frequency_neutral[[m]]
      ratio2 <-Size_of_beneficial_first_portion/mutation_frequency_neutral[[m]]
      ratio3 <-Size_of_beneficial_later_portion/mutation_frequency_neutral[[m]]
      ratio4 <-Size_of_with_known_beneficial_mutation/mutation_frequency_neutral[[m]]
      #print(paste("mutation ID ",m, "is neutral and portion of non-hitchhiker is ",ratio1,
      #           "      portion of beneficial first is ", ratio2,"    portion of beneficial later is ",ratio3))
      
      if(ratio1 < hitchhiker_cutoff){mutation_frequency_neutral_hitchhiker[[m]]=mutation_frequency_neutral[[m]];
      #within the hitchhikers, determine whether it comes before or after the beneficial mutation
      if(ratio2>ratio3){mutation_frequency_neutral_beneficial_first[[m]]=mutation_frequency_neutral[[m]];}
      else{ mutation_frequency_neutral_beneficial_later[[m]]=mutation_frequency_neutral[[m]];
      #print(paste("yes"))
      }
       if(ratio4/(1-ratio1)>known_cutoff){mutation_frequency_neutral_hitchhiker_with_known_ben[[m]]=mutation_frequency_neutral[[m]]}
      }
      else{mutation_frequency_neutral_not_hitchhiker[[m]]=mutation_frequency_neutral[[m]]}
      
      #mutation_frequency_neutral_hitchhiker[[m]]=mutation_frequency_neutral[[m]]
    }
    
  }
    if(mutation_frequency_neutral_hitchhiker[[m]]!=0|mutation_frequency_neutral_dev[[m]]!=0){
      
      if(mutation_frequency_neutral_hitchhiker[[m]]!=0){mutation_frequency_neutral_hitchhiker_developmental_union[[m]]=mutation_frequency_neutral_hitchhiker[[m]]}
      
      if(mutation_frequency_neutral_dev[[m]]!=0){mutation_frequency_neutral_hitchhiker_developmental_union[[m]]=mutation_frequency_neutral_dev[[m]]}
    }
  }
  
    
  
  
  
  if (End_Population_Size!=0){
    #turn all clonal sizes as frequency; particularly necessary when mean fitness update is not turned on: the entire population grows
    for (i in 1:length(mutation_frequency_beneficial)){
      mutation_frequency_beneficial[[i]]<-mutation_frequency_beneficial[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_beneficial_first_mutant)){
      mutation_frequency_beneficial_first_mutant[[i]]<-mutation_frequency_beneficial_first_mutant[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_beneficial_double_mutant)){
      mutation_frequency_beneficial_double_mutant[[i]]<-mutation_frequency_beneficial_double_mutant[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_neutral)){
      mutation_frequency_neutral[[i]]<-mutation_frequency_neutral[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_neutral_not_hitchhiker)){
      mutation_frequency_neutral_not_hitchhiker[[i]]<-mutation_frequency_neutral_not_hitchhiker[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_neutral_hitchhiker)){
      mutation_frequency_neutral_hitchhiker[[i]]<-mutation_frequency_neutral_hitchhiker[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_neutral_beneficial_first)){
      mutation_frequency_neutral_beneficial_first[[i]]<-mutation_frequency_neutral_beneficial_first[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_neutral_beneficial_later)){
      mutation_frequency_neutral_beneficial_later[[i]]<-mutation_frequency_neutral_beneficial_later[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_dev)){
      mutation_frequency_dev[[i]]<-mutation_frequency_dev[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_beneficial_dev)){
      mutation_frequency_beneficial_dev[[i]]<-mutation_frequency_beneficial_dev[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_neutral_dev)){
      mutation_frequency_neutral_dev[[i]]<-mutation_frequency_neutral_dev[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_neutral_hitchhiker_developmental_union)){
      mutation_frequency_neutral_hitchhiker_developmental_union[[i]]<-mutation_frequency_neutral_hitchhiker_developmental_union[[i]]/(2*End_Population_Size)
    }
    for (i in 1:length(mutation_frequency_neutral_hitchhiker_with_known_ben)){
      mutation_frequency_neutral_hitchhiker_with_known_ben[[i]]<-mutation_frequency_neutral_hitchhiker_with_known_ben[[i]]/(2*End_Population_Size)
    }
    
    
    mutation_frequency_beneficial<-mutation_frequency_beneficial[mutation_frequency_beneficial!=0]
    mutation_frequency_beneficial_first_mutant<-mutation_frequency_beneficial_first_mutant[mutation_frequency_beneficial_first_mutant!=0]
    mutation_frequency_beneficial_double_mutant<-mutation_frequency_beneficial_double_mutant[mutation_frequency_beneficial_double_mutant!=0]
    mutation_frequency_neutral<-mutation_frequency_neutral[mutation_frequency_neutral!=0]
    mutation_frequency_neutral_not_hitchhiker<-mutation_frequency_neutral_not_hitchhiker[mutation_frequency_neutral_not_hitchhiker!=0]
    mutation_frequency_neutral_hitchhiker<-mutation_frequency_neutral_hitchhiker[mutation_frequency_neutral_hitchhiker!=0]
    mutation_frequency_neutral_beneficial_first<-mutation_frequency_neutral_beneficial_first[mutation_frequency_neutral_beneficial_first!=0]
    mutation_frequency_neutral_beneficial_later<-mutation_frequency_neutral_beneficial_later[mutation_frequency_neutral_beneficial_later!=0]
    mutation_frequency_dev <-mutation_frequency_dev[mutation_frequency_dev!=0]
    mutation_frequency_beneficial_dev <-mutation_frequency_beneficial_dev[mutation_frequency_beneficial_dev!=0]
    mutation_frequency_neutral_dev <-mutation_frequency_neutral_dev[mutation_frequency_neutral_dev!=0]
    mutation_frequency_neutral_hitchhiker_developmental_union <- mutation_frequency_neutral_hitchhiker_developmental_union[mutation_frequency_neutral_hitchhiker_developmental_union!=0]
    mutation_frequency_neutral_hitchhiker_with_known_ben<-mutation_frequency_neutral_hitchhiker_with_known_ben[mutation_frequency_neutral_hitchhiker_with_known_ben!=0]
    
    
  }
  
  
  # print(paste('the ratio of all hitchhikers to those hitchhiking with known ben. mutations is',
  #             length(mutation_frequency_neutral_hitchhiker)/length(mutation_frequency_neutral_hitchhiker_with_known_ben)))
  # print(paste('the ratio of all beneficial mutations to known ben. mutations is',
  #             length(beneficial_mut_ID)/length(known_beneficial_mut_ID)))
  
  Total_mutation_frequency_beneficial<-list.append(Total_mutation_frequency_beneficial,mutation_frequency_beneficial)
  Total_mutation_frequency_beneficial_first_mutant<-list.append(Total_mutation_frequency_beneficial_first_mutant,mutation_frequency_beneficial_first_mutant)
  Total_mutation_frequency_beneficial_double_mutant<-list.append(Total_mutation_frequency_beneficial_double_mutant,mutation_frequency_beneficial_double_mutant)
  Total_mutation_frequency_neutral<-list.append(Total_mutation_frequency_neutral,mutation_frequency_neutral)
  Total_mutation_frequency_neutral_not_hitchhiker<-list.append(Total_mutation_frequency_neutral_not_hitchhiker,mutation_frequency_neutral_not_hitchhiker)
  Total_mutation_frequency_neutral_hitchhiker<-list.append(Total_mutation_frequency_neutral_hitchhiker,mutation_frequency_neutral_hitchhiker)
  Total_mutation_frequency_neutral_beneficial_first<-list.append(Total_mutation_frequency_neutral_beneficial_first,mutation_frequency_neutral_beneficial_first)
  Total_mutation_frequency_neutral_beneficial_later<-list.append(Total_mutation_frequency_neutral_beneficial_later,mutation_frequency_neutral_beneficial_later)
  Total_mutation_frequency_dev <-list.append(Total_mutation_frequency_dev, mutation_frequency_dev)
  Total_mutation_frequency_beneficial_dev <-list.append(Total_mutation_frequency_beneficial_dev, mutation_frequency_beneficial_dev)
  Total_mutation_frequency_neutral_dev <-list.append(Total_mutation_frequency_neutral_dev, mutation_frequency_neutral_dev)
  Total_mutation_frequency_neutral_hitchhiker_developmental_union <-list.append(Total_mutation_frequency_neutral_hitchhiker_developmental_union, mutation_frequency_neutral_hitchhiker_developmental_union)
  Total_mutation_frequency_neutral_hitchhiker_with_known_ben<-list.append(Total_mutation_frequency_neutral_hitchhiker_with_known_ben, mutation_frequency_neutral_hitchhiker_with_known_ben)
  print(paste("This is patient number ",N, " has number of mutations", last_mutation_ID ))
  print(paste("The ratio of end population to starting population is ", End_Population_Size/Population_size, 'and end of development occurred at', End_of_development, 'after mutation', last_development_mutation_ID ))
  
  
  
}                                                  #end of number of participant loop


end_time <- Sys.time()

code_ran_for<- end_time - start_time

code_ran_for<-round(code_ran_for,1)


##################################    OVERVIEW of ALL PATIENTS    ###################################################
print(paste(' number of non-extinct lineages is', num_non_extinct_lineage))
# 
Total_mutation_frequency_beneficial<-unlist(Total_mutation_frequency_beneficial)
Total_mutation_frequency_beneficial_first_mutant<-unlist(Total_mutation_frequency_beneficial_first_mutant)
Total_mutation_frequency_beneficial_double_mutant<-unlist(Total_mutation_frequency_beneficial_double_mutant)
Total_mutation_frequency_neutral<-unlist(Total_mutation_frequency_neutral)
Total_mutation_frequency_neutral_not_hitchhiker<-unlist(Total_mutation_frequency_neutral_not_hitchhiker)
Total_mutation_frequency_neutral_hitchhiker<-unlist(Total_mutation_frequency_neutral_hitchhiker)
Total_mutation_frequency_neutral_beneficial_first<-unlist(Total_mutation_frequency_neutral_beneficial_first)
Total_mutation_frequency_neutral_beneficial_later<-unlist(Total_mutation_frequency_neutral_beneficial_later)
Total_mutation_frequency_dev <- unlist(Total_mutation_frequency_dev)
Total_mutation_frequency_beneficial_dev <- unlist(Total_mutation_frequency_beneficial_dev)
Total_mutation_frequency_neutral_dev <- unlist(Total_mutation_frequency_neutral_dev)
Total_mutation_frequency_neutral_hitchhiker_developmental_union<-unlist(Total_mutation_frequency_neutral_hitchhiker_developmental_union)
Total_mutation_frequency_neutral_hitchhiker_with_known_ben<-unlist(Total_mutation_frequency_neutral_hitchhiker_with_known_ben)
# 
# Total_mutation_frequency_beneficial<-Total_mutation_frequency_beneficial[Total_mutation_frequency_beneficial!=0]
# Total_mutation_frequency_beneficial_first_mutant<-Total_mutation_frequency_beneficial_first_mutant[Total_mutation_frequency_beneficial_first_mutant!=0]
# Total_mutation_frequency_beneficial_double_mutant<-Total_mutation_frequency_beneficial_double_mutant[Total_mutation_frequency_beneficial_double_mutant!=0]
# Total_mutation_frequency_neutral<-Total_mutation_frequency_neutral[Total_mutation_frequency_neutral!=0]
# Total_mutation_frequency_neutral_not_hitchhiker<-Total_mutation_frequency_neutral_not_hitchhiker[Total_mutation_frequency_neutral_not_hitchhiker!=0]
# Total_mutation_frequency_neutral_hitchhiker<-Total_mutation_frequency_neutral_hitchhiker[Total_mutation_frequency_neutral_hitchhiker!=0]
# Total_mutation_frequency_neutral_beneficial_first<-Total_mutation_frequency_neutral_beneficial_first[Total_mutation_frequency_neutral_beneficial_first!=0]
# Total_mutation_frequency_neutral_beneficial_later<-Total_mutation_frequency_neutral_beneficial_later[Total_mutation_frequency_neutral_beneficial_later!=0]
# Total_mutation_frequency_dev <-Total_mutation_frequency_dev[Total_mutation_frequency_dev!=0]
# Total_mutation_frequency_beneficial_dev <-Total_mutation_frequency_beneficial_dev[Total_mutation_frequency_beneficial_dev!=0]
# Total_mutation_frequency_neutral_dev <-Total_mutation_frequency_neutral_dev[Total_mutation_frequency_neutral_dev!=0]
# Total_mutation_frequency_neutral_hitchhiker_developmental_union <- Total_mutation_frequency_neutral_hitchhiker_developmental_union[Total_mutation_frequency_neutral_hitchhiker_developmental_union!=0]
# Total_mutation_frequency_neutral_hitchhiker_with_known_ben<-Total_mutation_frequency_neutral_hitchhiker_with_known_ben[Total_mutation_frequency_neutral_hitchhiker_with_known_ben!=0]

if (length(Total_mutation_frequency_neutral_hitchhiker_developmental_union)==
    length(Total_mutation_frequency_neutral_hitchhiker)+length(Total_mutation_frequency_neutral_dev))
{print(paste('the union is a union'))}

#save the data 

write.csv(Total_mutation_frequency_beneficial, file = "beneficial_mutation_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_beneficial_first_mutant, file = "beneficial_mutation_first_mutant_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_beneficial_double_mutant, file = "beneficial_mutation_double_mutant_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_neutral, file = "neutral_mutation_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_neutral_not_hitchhiker, file = "neutral_mutation_not_hitchhiker_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_neutral_beneficial_first, file = "neutral_mutation_beneficial_first_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_neutral_beneficial_later, file = "neutral_mutation_hitchhiker_beneficial_later_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_neutral_hitchhiker, file = "neutral_mutation_hitchhiker_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_dev, file = "dev_mutation_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_beneficial_dev, file = "beneficial_dev_mutation_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_neutral_dev, file = "neutral_dev_mutation_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
write.csv(Total_mutation_frequency_neutral_hitchhiker_developmental_union, file = "neutral_dev_hitchhiker_union_r_1p2_lower_dev_rate_uben_single_s_10p_70_gen_run_11.csv")
# write.csv(Total_mutation_frequency_neutral_hitchhiker_with_known_ben, file = "M1_neutral_hitchhiker_with_known_ben_r_120p_s_DFE_smax_14p_percent_70_gen_run_8.csv")
