# 21/11/2018
# k-hit
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

#par(mfrow=c(3,1))

# Bio_system <-list(list(list(list(1,0)),10^5));#clone_info<- (mutation, clone size)
#starting from one single clone with three mutations
#Bio_system <-list(list(list(list(1,0.01),list(2,0.02),list(3,0.0015)),100),list(list(list(4,0.1)),100000));#clone_info<- (mutation, clone size)
    #starting from one single clone with three mutations

#Define function 'Divide'
Divide<- function(clone_info, dt, B0, D0,mean_fitness){
  #dt = a small interval of time, 
  #D0 = symmetric rate to differentiated
  #B0 = symmetric rate to self renewal
  #B  = modified birth rate including fitness advantage
  
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


#Define function 'mutation_fitness' 
mutation_fitness <- function(DFE,ratio_neutral_beneficial,s1) {
  
  random_number<- runif(1, 0, 1)
  threshold= 1/(1+ratio_neutral_beneficial)
  if(random_number<threshold){
    
    s=s1
  
  }else{s=0}
  
  return(s)                                  
}


#Define function 'mutate'
mutate<-function(clone_info, dt, u, last_mutation_ID,ratio_neutral_beneficial,s1){ #u=mutation rate
  mutations=clone_info[[1]] #[[ID1, s1],[ID2,s2]]
  clone_size=clone_info[[2]]
  
  # if(clone_size*u*dt>1){                                       #do not mutate clones with size<1/U
  
  
  number_of_mutations <-rpois(1, clone_size*u*dt) ;
     #print(paste("mean of poisson =", clone_size*u*dt))
      # print(paste("number of mutations =", number_of_mutations))

list_of_new_clones=list()
if(number_of_mutations!=0){
for (i in 1:number_of_mutations){
  last_mutation_ID = last_mutation_ID+1
        #print(paste("last mutation ID =", last_mutation_ID))
  new_fitness_effect = mutation_fitness(1,ratio_neutral_beneficial,s1)
#         print("new fitness effect =", new_fitness_effect)
  new_mutation<-list(last_mutation_ID, new_fitness_effect)
  new_clone_info <-list(list.append(mutations,new_mutation),1)    #[[[ID1,s1],[ID2,s2],....],n]
#         print("new clone =", new_clone)
  list_of_new_clones<-list.append(list_of_new_clones,new_clone_info)
#         print("list of new clones =", list_of_new_clones)
}
}   
   
   
   
  # }else{
  #   
  #   list_of_new_clones<-list()
  #   last_mutation_ID = last_mutation_ID
  #   number_of_mutations=0
  #   
  # }
  
result<-list(list_of_new_clones, last_mutation_ID, number_of_mutations)     
return (result)                                      #result[1] returns the list of new clones
}
#what happens when number_of_mutations=0?

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
Population_size=10^5                                 #is the starting population size of clone
hitchhiker_cutoff=0.5;
# benefit=0.04;
dt=0.1;#in units of generation
list_of_ages = list(21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 37.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 45.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5,
                    69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5, 73.5)

Number_of_Participants = length(list_of_ages)
Total_mutation_frequency_beneficial<-list()
Total_mutation_frequency_beneficial_first_mutant<-list()
Total_mutation_frequency_beneficial_double_mutant<-list()
Total_mutation_frequency_neutral<-list()
Total_mutation_frequency_neutral_not_hitchhiker<-list()
Total_mutation_frequency_neutral_hitchhiker<-list()
Total_mutation_frequency_neutral_beneficial_first<-list()
Total_mutation_frequency_neutral_beneficial_later<-list()
Clonal_histories_per_person<-list()
Entire_experiment<-list()

cat(paste("largest_driver_vaf","        ","largest_neutral_vaf","        ", "cooccurrence"),file="simulated_correlation_plot_Martincorena_ages_unseen_drivers.txt",sep="\n")

unseen_fold = 1
denom = unseen_fold+1
cooccurred_cd_top_driver_larger_than_top_neutral=0
cooccurred_cd_top_driver_smaller_than_top_neutral=0
times_top_driver_larger_than_top_neutral=0
times_top_driver_smaller_than_top_neutral=0
top_driver_larger_than_top_neutral_cd_cooccurred=0
top_driver_smaller_than_top_neutral_cd_cooccurred=0
top_driver_larger_than_top_neutral_cd_not_cooccurred=0
top_driver_smaller_than_top_neutral_cd_not_cooccurred=0
times_cooccurred=0
times_not_cooccurred=0
top_driver_vaf<-list()
top_neutral_vaf<-list()
for (N in 1:Number_of_Participants) {
  
number_of_generations=list_of_ages[[N]];#number of generations
lifespan=number_of_generations/dt; #number of runs
Bio_system <-list(list(list(list(1,0)),Population_size));
#population_limit=10^15;
t=0;
last_mutation_ID=1;
B0=1;
D0=1;
s1=0.11
u1=3.3*10^{-5}
u2=1.2*10^{-5}
u_ben=(u1+u2)*(unseen_fold+1)
u_neu=2.9*10^-4
u = u_ben + u_neu                      #mutation rate (neutral and beneficial)
ratio_neutral_beneficial=u_neu/u_ben;           #R=1 for beneficial mutation
mean_fitness=0;
End_Population_Size=0;
Clonal_histories<-list() 

Stored_mean_fitness<-list()
Stored_number_of_clones<-list()
Bio_system_record<-list()
mean_fitness_record<-list()
Number_of_nonextinct_clones_record<-list()

for (t in 1:lifespan) {


#mutate Bio_system 
  result=list();
  list_of_new_clones=list()
  number_of_mutations=list()

          for (m in 1:length(Bio_system)){
            
            if(Bio_system[[m]][[2]]!=0){                            #only mutate non-extinct clones

          result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,s1) #consider the first clone
          list_of_new_clones[[m]]<- result[[m]][[1]]
          last_mutation_ID <- result[[m]][[2]]
          number_of_mutations[[m]]<- result[[m]][[3]]
             #only call the function 'mutate' once for each clone
          
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
            
              if(TRUE){
            Bio_system[[i]][[2]]=Divide(Bio_system[[i]], dt, B0, D0,mean_fitness)[[1]]   #update clone size
                   }    #else{Bio_system[[i]][[2]]=population_limit}
            
            }                                                       #consider relative fitness
            Clonal_sizes[[i]]=Bio_system[[i]][[2]];
          
          }
  #print("Finished dividing")
          
#update mean fitness
          result2<-update_mean_fitness(Bio_system)
          mean_fitness<-result2[[1]]
          Stored_mean_fitness<-list.append(Stored_mean_fitness,mean_fitness)
          
          total_pop<-result2[[2]]


          
          t=t+dt;
          #print(paste("this is generation ",t))
         
}






#######################################START OF PART FOR CLONAL HISTORY#################################
# #to make all vectors of Clonal_histories the same length
# T<- lifespan-1
# 
# for (t in 1:T) {
# 
#   List_of_Clonal_histories<-Clonal_histories[[t]]
#   Number_of_zeros_to_add <-length(Bio_system)-length(List_of_Clonal_histories)
#   for (j in 1:Number_of_zeros_to_add) {
# 
# 
#     List_of_Clonal_histories<-list.append(List_of_Clonal_histories,0)
# 
#   }
#   Clonal_histories[[t]]<-List_of_Clonal_histories
# }
# 
# 
# 
# #make data frame of clonal histories
# 
# Clone_history<-list()
# Total_history<-c()
#         #change from longitudinal to transverse
#         for (i in 1:length(Bio_system)) {
#         
#            for (t in 1:lifespan){                                        
#         
#             Clone_history[[t]]<- Clonal_histories[[t]][[i]]   #is the clonal history of the i-th clone
#         
#            }
#             cln<-unlist(Clone_history)
#             Total_history <- cbind(Total_history,cln)
#             #Total_history[[a,b]] is the clonal size of the b-th clone in generation a
# 
#           }
# 
# Clonal_histories_per_person[[N]]<-Total_history
#stores clonal histories of each person into a separate entry



# mean_fitness_record[[N]]<-Stored_mean_fitness
# #stores the the mean fitness trajectory of each patient into a separate entry
# Number_of_nonextinct_clones_record[[N]]<-Stored_number_of_clones
# #stores the record of the number of nonextinct clones of each patient into a separate entry
# Entire_experiment[[N]]<-Bio_system_record
# #stores the Bio_system_record of each patient into a separate entry


# 
# #implement this only for efficiency!!!!
# if(FALSE){
# Total_history_efficiency<- Total_history
# useless_intermediate<-list()
# for (b in 1:length(Bio_system)) {
#   
# 
# if(Total_history_efficiency[[lifespan,b]]!=0){
#   useless_intermediate<-cbind.(useless_intermediate,Total_history_efficiency[[b]])
#                                               }
#   
#                                 }
# Total_history_efficiency <- useless_intermediate  
# }
#is the total history of non-extinct clones
#######################################END OF PART FOR CLONAL HISTORY##############################################

#find mutation frequency for this patient
mutation_frequency_beneficial<-list()                              #all beneficial mutations
mutation_frequency_beneficial_first_mutant<-list()
mutation_frequency_beneficial_double_mutant<-list()
mutation_frequency_neutral<-list()                                #all neutral mutations
mutation_frequency_neutral_not_hitchhiker<-list()   
mutation_frequency_neutral_hitchhiker<-list() 
mutation_frequency_neutral_beneficial_first<-list()
mutation_frequency_neutral_beneficial_later<-list()
this_mutation_cooccurs_with_drivers<-list()
this_hitchhiker_cooccurs_with_drivers<-list()

for (n in 1:length(Bio_system)) {

  Size<- Bio_system[[n]][[2]]
  if(Size!=0){End_Population_Size=End_Population_Size+Size}
  
}

for (m in 1:last_mutation_ID){mutation_frequency_beneficial[[m]]=0;
                              mutation_frequency_beneficial_first_mutant[[m]]=0;
                              mutation_frequency_beneficial_double_mutant[[m]]=0;
                              mutation_frequency_neutral[[m]]=0;
                              mutation_frequency_neutral_not_hitchhiker[[m]]=0;
                              mutation_frequency_neutral_hitchhiker[[m]]=0;
                              mutation_frequency_neutral_beneficial_first[[m]]=0;
                              mutation_frequency_neutral_beneficial_later[[m]]=0;
                              this_mutation_cooccurs_with_drivers[[m]]=0
                              this_hitchhiker_cooccurs_with_drivers[[m]]=0
                              }    #initialization

for (m in 2:last_mutation_ID) {
  Size_of_non_hitchhiker_portion=0;
  Size_of_beneficial_first_portion=0;
  Size_of_beneficial_later_portion=0;
  this_mutation_cooccurs_with_drivers[[m]]<-list()
  
  for (n in 1:length(Bio_system)) {
    
    CloneID<- Bio_system[[n]][[1]]
    Size<- Bio_system[[n]][[2]]
    
    
    
    Type_of_neutral_mutationID=0;                           #this can be "Not Hitchhiker", "Beneficial first" or "Beneficial later" 
    where_first_beneficial_mutation_comes_in=0;           #to be the position of the 1st beneficial mutation in the n-th clone
    number_of_beneficial_mutations=0;
    total_fitness=0;
    
    #calculate total fitness of clone
    for (k in 1:length(CloneID)) {
      total_fitness=total_fitness+CloneID[[k]][[2]];
    }
    
    #determine if the clone is neutral
    if(total_fitness==0){                                #relies on the fact that WT mutation is neutral
      Type_of_neutral_mutationID="Not Hitchhiker"
       #print(paste("not hitchhiker"))
      
    }else{                     
      #if clone not neutral, then find when the first beneficial mutation comes in
      for (k in 1:length(CloneID)) {
        if(CloneID[[k]][[2]]!=0){
      where_first_beneficial_mutation_comes_in=k;         #position of the 1st beneficial mutation   
            #print(paste("where_first_beneficial_mutation_comes_in is",where_first_beneficial_mutation_comes_in))
      break
                                }
      
      }
      
      for (k in 1:length(CloneID)){
        if (CloneID[[k]][[2]]!=0){
          this_mutation_cooccurs_with_drivers[[m]]<- list.append(this_mutation_cooccurs_with_drivers[[m]], CloneID[[k]][[1]])
        }
      }
      
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
        
        if(CloneID[[k]][[2]]!=0){
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
          this_hitchhiker_cooccurs_with_drivers[[m]] <- this_mutation_cooccurs_with_drivers[[m]]
          
          if(k<where_first_beneficial_mutation_comes_in & where_first_beneficial_mutation_comes_in!=0){ Type_of_neutral_mutationID="Beneficial later"}
          if(k>where_first_beneficial_mutation_comes_in & where_first_beneficial_mutation_comes_in!=0){ Type_of_neutral_mutationID="Beneficial first"}
          if(where_first_beneficial_mutation_comes_in==0){ Type_of_neutral_mutationID="Not Hitchhiker"}
                       #print(paste("The neutral mutation ", m," belongs to the type", Type_of_neutral_mutationID, "in clone ",n))
          
          
          
          if( Type_of_neutral_mutationID=="Not Hitchhiker"){
            Size_of_non_hitchhiker_portion=Size_of_non_hitchhiker_portion+Size;}
          if( Type_of_neutral_mutationID=="Beneficial first"){
            Size_of_beneficial_first_portion=Size_of_beneficial_first_portion+Size;}
          if( Type_of_neutral_mutationID=="Beneficial later"){
            Size_of_beneficial_later_portion=Size_of_beneficial_later_portion+Size} 
              }
        
           break                  }        #no need to look at the rest of the CloneID since a mutation only occurs once
      
      
      
    }
    
  }
  if(mutation_frequency_neutral[[m]]!=0){  #it is possible that that a mutation is found in two clones, a neutral clone and a hitchhiker clone
  ratio1 <-Size_of_non_hitchhiker_portion/mutation_frequency_neutral[[m]]
  ratio2 <-Size_of_beneficial_first_portion/mutation_frequency_neutral[[m]]
  ratio3 <-Size_of_beneficial_later_portion/mutation_frequency_neutral[[m]]
     #print(paste("mutation ID ",m, "is neutral and portion of non-hitchhiker is ",ratio1,
      #           "      portion of beneficial first is ", ratio2,"    portion of beneficial later is ",ratio3))
  
  if(ratio1 < hitchhiker_cutoff){mutation_frequency_neutral_hitchhiker[[m]]=mutation_frequency_neutral[[m]];
                                #within the hitchhikers, determine whether it comes before or after the beneficial mutation
                                if(ratio2>ratio3){mutation_frequency_neutral_beneficial_first[[m]]=mutation_frequency_neutral[[m]];}
                                                 else{ mutation_frequency_neutral_beneficial_later[[m]]=mutation_frequency_neutral[[m]];
                                                      #print(paste("yes"))
                                                     }

                                }
  else{mutation_frequency_neutral_not_hitchhiker[[m]]=mutation_frequency_neutral[[m]]}
  
  #mutation_frequency_neutral_hitchhiker[[m]]=mutation_frequency_neutral[[m]]
  }
  # if(mutation_frequency[[m]]!=0){
  #   #print(paste("The mutation frequency is ",mutation_frequency[[m]],"for mutation ID ",m))
  # }
}

#turn all clonal sizes as frequency; particularly necessary when mean fitness update is not turned on: the entire population grows
for (i in 1:length(mutation_frequency_beneficial)){
mutation_frequency_beneficial[[i]]<-mutation_frequency_beneficial[[i]]/End_Population_Size
}
for (i in 1:length(mutation_frequency_beneficial_first_mutant)){
  mutation_frequency_beneficial_first_mutant[[i]]<-mutation_frequency_beneficial_first_mutant[[i]]/End_Population_Size
}
for (i in 1:length(mutation_frequency_beneficial_double_mutant)){
  mutation_frequency_beneficial_double_mutant[[i]]<-mutation_frequency_beneficial_double_mutant[[i]]/End_Population_Size
}
for (i in 1:length(mutation_frequency_neutral)){
  mutation_frequency_neutral[[i]]<-mutation_frequency_neutral[[i]]/End_Population_Size
}
for (i in 1:length(mutation_frequency_neutral_not_hitchhiker)){
  mutation_frequency_neutral_not_hitchhiker[[i]]<-mutation_frequency_neutral_not_hitchhiker[[i]]/End_Population_Size
}
for (i in 1:length(mutation_frequency_neutral_hitchhiker)){
  mutation_frequency_neutral_hitchhiker[[i]]<-mutation_frequency_neutral_hitchhiker[[i]]/End_Population_Size
}
for (i in 1:length(mutation_frequency_neutral_beneficial_first)){
  mutation_frequency_neutral_beneficial_first[[i]]<-mutation_frequency_neutral_beneficial_first[[i]]/End_Population_Size
}
for (i in 1:length(mutation_frequency_neutral_beneficial_later)){
  mutation_frequency_neutral_beneficial_later[[i]]<-mutation_frequency_neutral_beneficial_later[[i]]/End_Population_Size
}

largest_neutral_vaf<-max(unlist(mutation_frequency_neutral))
largest_neutral_mutation_ID <-which.max(mutation_frequency_neutral)
largest_neutral_hitchhikes_with_drivers_whose_ID <- this_hitchhiker_cooccurs_with_drivers[[largest_neutral_mutation_ID]]
largest_driver_vaf<-max(unlist(mutation_frequency_beneficial))
largest_driver_mutation_ID <-which.max(mutation_frequency_beneficial)

tmp<-list()
for (number in 1:length(mutation_frequency_beneficial)){
  
  if (number %% denom ==0){
    tmp <- list.append(tmp, mutation_frequency_beneficial[[number]])
  }
}

largest_seen_driver_vaf <- max(unlist(tmp))
largest_seen_driver_mutation_ID <-which.max(tmp)*denom

if (largest_neutral_vaf<largest_driver_vaf){
  times_top_driver_larger_than_top_neutral=times_top_driver_larger_than_top_neutral+1
  if (is.element(largest_driver_mutation_ID, largest_neutral_hitchhikes_with_drivers_whose_ID)==TRUE){
    cooccurred_cd_top_driver_larger_than_top_neutral=cooccurred_cd_top_driver_larger_than_top_neutral+1
  }
}


if (largest_neutral_vaf>largest_driver_vaf){
  times_top_driver_smaller_than_top_neutral=times_top_driver_smaller_than_top_neutral+1
  if (is.element(largest_driver_mutation_ID, largest_neutral_hitchhikes_with_drivers_whose_ID)==TRUE){
    cooccurred_cd_top_driver_smaller_than_top_neutral=cooccurred_cd_top_driver_smaller_than_top_neutral+1
  }
}
if (is.element(largest_driver_mutation_ID, largest_neutral_hitchhikes_with_drivers_whose_ID)==TRUE){
  cooccurrence = 'True'
  times_cooccurred = times_cooccurred+1
  if (largest_neutral_vaf<largest_driver_vaf){
    top_driver_larger_than_top_neutral_cd_cooccurred=top_driver_larger_than_top_neutral_cd_cooccurred+1}
  if (largest_neutral_vaf>largest_driver_vaf){
    top_driver_smaller_than_top_neutral_cd_cooccurred=top_driver_smaller_than_top_neutral_cd_cooccurred+1
  }
}else{
  cooccurrence = 'False'
  times_not_cooccurred = times_not_cooccurred + 1
  if (largest_neutral_vaf<largest_driver_vaf){
    top_driver_larger_than_top_neutral_cd_not_cooccurred=top_driver_larger_than_top_neutral_cd_not_cooccurred+1}
  if (largest_neutral_vaf>largest_driver_vaf){
    top_driver_smaller_than_top_neutral_cd_not_cooccurred=top_driver_smaller_than_top_neutral_cd_not_cooccurred+1
  }
}
top_driver_vaf<-list.append(top_driver_vaf,largest_driver_vaf)
top_neutral_vaf<-list.append(top_neutral_vaf,largest_neutral_vaf)

if (is.element(largest_seen_driver_mutation_ID, largest_neutral_hitchhikes_with_drivers_whose_ID)==TRUE){
  cooccurrence_with_seen_drivers = 'True'}else{cooccurrence_with_seen_drivers = 'False'}
cat(paste(largest_seen_driver_vaf, "        ",largest_neutral_vaf,"        ", cooccurrence_with_seen_drivers),file="simulated_correlation_plot_Martincorena_ages_unseen_drivers_fold_1.txt",append=TRUE, sep="\n")
cooccurrence = 'error'

# Total_mutation_frequency_beneficial<-list.append(Total_mutation_frequency_beneficial,mutation_frequency_beneficial)
# Total_mutation_frequency_beneficial_first_mutant<-list.append(Total_mutation_frequency_beneficial_first_mutant,mutation_frequency_beneficial_first_mutant)
# Total_mutation_frequency_beneficial_double_mutant<-list.append(Total_mutation_frequency_beneficial_double_mutant,mutation_frequency_beneficial_double_mutant)
# Total_mutation_frequency_neutral<-list.append(Total_mutation_frequency_neutral,mutation_frequency_neutral)
# Total_mutation_frequency_neutral_not_hitchhiker<-list.append(Total_mutation_frequency_neutral_not_hitchhiker,mutation_frequency_neutral_not_hitchhiker)
# Total_mutation_frequency_neutral_hitchhiker<-list.append(Total_mutation_frequency_neutral_hitchhiker,mutation_frequency_neutral_hitchhiker)
# Total_mutation_frequency_neutral_beneficial_first<-list.append(Total_mutation_frequency_neutral_beneficial_first,mutation_frequency_neutral_beneficial_first)
# Total_mutation_frequency_neutral_beneficial_later<-list.append(Total_mutation_frequency_neutral_beneficial_later,mutation_frequency_neutral_beneficial_later)
                   print(paste("This is patient number ",N, " has number of mutations", last_mutation_ID ))
                   print(paste("The ratio of end population to starting population is ", End_Population_Size/Population_size ))
                   

}                                                  #end of number of participant loop

print(paste('likelihood of co-occurrence conditioned
on observing top driver vaf larger than top neutral vaf', cooccurred_cd_top_driver_larger_than_top_neutral/times_top_driver_larger_than_top_neutral))
print(paste('likelihood of co-occurrence conditioned
            on observing top driver vaf smaller than top neutral vaf', cooccurred_cd_top_driver_smaller_than_top_neutral/times_top_driver_smaller_than_top_neutral))
print(paste('likelihood of observing top driver vaf larger than top neutral vaf conditioned
            on cooccurrence', top_driver_larger_than_top_neutral_cd_cooccurred/times_cooccurred))
print(paste('likelihood of observing top driver vaf smaller than top neutral vaf conditioned
            on cooccurrence', top_driver_smaller_than_top_neutral_cd_cooccurred/times_cooccurred))
print(paste('likelihood of observing top driver vaf larger than top neutral vaf conditioned
            on cooccurrence', top_driver_larger_than_top_neutral_cd_not_cooccurred/times_not_cooccurred))
print(paste('likelihood of observing top driver vaf smaller than top neutral vaf conditioned
            on cooccurrence', top_driver_smaller_than_top_neutral_cd_not_cooccurred/times_not_cooccurred))

end_time <- Sys.time()

code_ran_for<- end_time - start_time

code_ran_for<-round(code_ran_for,1)

######################################################################################################################
###################################     PLOT HISTORY OF A SINGLE PATIENT     ###################################


# 
# M=1     #specify which patient's history
# 
# #plot the clonal histories
# 
# History<-Clonal_histories_per_person[[M]]
# timescale<-c(1:lifespan)
# matplot(timescale, log(History[,1:length(Bio_system)]),col=1:length(Bio_system),lty=2, lwd=1,
#         type = "l",
#         xlab="number of generations",
#         ylab="Clone Size",
#         xlim=c(0,lifespan),
#        #ylim=c(0,110000)
# 
# )
# #title(main = "Trajectories of mutants, 1000 generations, 8000 runs")
# 
# 
# 
# #plot mean fitness with time
# 
# Record_mean_fitness<-unlist(mean_fitness_record[[M]])
# 
# plot(timescale,Record_mean_fitness,xlab="number of generations",
#      ylab="Mean fitness",lty=2)
# 
# 
# 
# #plot mean fitness with time
# 
# Number_of_nonextinct_clones<-unlist(Number_of_nonextinct_clones_record[[M]])
# 
# plot(timescale,Number_of_nonextinct_clones,xlab="number of generations",
#      ylab="Number of nonextinct Clones",lty=2)
# 
# 
# 
# 
# grid.arrange(plot1,plot2,plot3, nrow = 3)
# 



##################################    OVERVIEW of ALL PATIENTS    ###################################################
Total_mutation_frequency_beneficial<-unlist(Total_mutation_frequency_beneficial)
Total_mutation_frequency_beneficial_first_mutant<-unlist(Total_mutation_frequency_beneficial_first_mutant)
Total_mutation_frequency_beneficial_double_mutant<-unlist(Total_mutation_frequency_beneficial_double_mutant)
Total_mutation_frequency_neutral<-unlist(Total_mutation_frequency_neutral)
Total_mutation_frequency_neutral_not_hitchhiker<-unlist(Total_mutation_frequency_neutral_not_hitchhiker)
Total_mutation_frequency_neutral_hitchhiker<-unlist(Total_mutation_frequency_neutral_hitchhiker)
Total_mutation_frequency_neutral_beneficial_first<-unlist(Total_mutation_frequency_neutral_beneficial_first)
Total_mutation_frequency_neutral_beneficial_later<-unlist(Total_mutation_frequency_neutral_beneficial_later)

Total_mutation_frequency_beneficial<-Total_mutation_frequency_beneficial[Total_mutation_frequency_beneficial!=0]
Total_mutation_frequency_beneficial_first_mutant<-Total_mutation_frequency_beneficial_first_mutant[Total_mutation_frequency_beneficial_first_mutant!=0]
Total_mutation_frequency_beneficial_double_mutant<-Total_mutation_frequency_beneficial_double_mutant[Total_mutation_frequency_beneficial_double_mutant!=0]
Total_mutation_frequency_neutral<-Total_mutation_frequency_neutral[Total_mutation_frequency_neutral!=0]
Total_mutation_frequency_neutral_not_hitchhiker<-Total_mutation_frequency_neutral_not_hitchhiker[Total_mutation_frequency_neutral_not_hitchhiker!=0]
Total_mutation_frequency_neutral_hitchhiker<-Total_mutation_frequency_neutral_hitchhiker[Total_mutation_frequency_neutral_hitchhiker!=0]
Total_mutation_frequency_neutral_beneficial_first<-Total_mutation_frequency_neutral_beneficial_first[Total_mutation_frequency_neutral_beneficial_first!=0]
Total_mutation_frequency_neutral_beneficial_later<-Total_mutation_frequency_neutral_beneficial_later[Total_mutation_frequency_neutral_beneficial_later!=0]







