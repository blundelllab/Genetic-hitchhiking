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
mutation_fitness <- function(DFE,ratio_neutral_beneficial,benefit) {
  
  random_number<- runif(1, 0, 1)
  threshold= 1/(1+ratio_neutral_beneficial)
  if(random_number<threshold){
    
    b = 2
    d = 0.06
    smax = 1
    smin = 0
    
    DFE <- function(x) exp(-(x**2/d**2)**(b/2))
    
    dist_DFE <- distr::AbscontDistribution(d=DFE)
    fitness <- distr::r(dist_DFE)
    
    repeat{
      s=fitness(1)
      if(s < smax & s > smin){
        break
      }
      
    }
    
  }else{s=0}
  
  return(s)                                  
}


#Define function 'mutate'
mutate<-function(clone_info, dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit){ #u=mutation rate
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
  new_fitness_effect = mutation_fitness(1,ratio_neutral_beneficial,benefit)
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
list_of_ages_all = list(67.38946, 63.04175, 58.91307, 53.31417, 40.29295, 54.52977, 66.940453, 
                        71.600273, 83.09925, 64.10677, 65.83984, 35.10199, 58.17659, 67.24983, 72.07666, 59.274471, 48.67899, 66.25599, 85.563316, 65.097878, 
                        64.62423, 55.26078, 82.724159, 74.75428, 67.12936, 65.49487, 74.03696, 80.950035, 70.77071, 77.563316, 80.75839, 73.71116, 81.472961, 75.28542, 61.18823, 64.495552, 72.04107, 61.316906, 50.11088, 77.21835, 67.559204, 54.64476, 76.53388, 69.82341, 66.160164, 53.05133, 68.427101, 54.10541, 85.96577, 55.827515, 76.695412, 58.44764, 62.959618, 76.793976, 76.84326, 77.801506, 81.300476, 24.82136, 79.27721, 82.603699, 78.346336, 77.221085, 60.944557, 66.362762, 59.64134, 67.57016, 79.430527, 50.132786, 64.982887, 67.61396, 56.49829, 66.234085, 61.2731, 71.088295, 60.804928, 57.95483, 73.92744, 69.059547, 56.095825, 65.505821, 73.06503, 65.83162, 54.91034, 52.818619, 68.89254, 81.15263, 65.51403, 55.24983, 69.55236, 75.843941, 75.167694, 77.31964, 53.7577, 74.54073, 74.54346, 67.077347, 23.83299, 82.817245, 76.14236, 72.33675, 66.989731, 78.85558, 39.90965, 67.28542, 73.735794, 78.266937, 73.492126, 83.18412, 60.175224, 69.3306, 64.342232, 65.185486, 71.84942, 62.66393, 65.730324, 73.36345, 57.03765, 48.40794, 69.503082, 74.29158, 68.67078, 65.24298, 61.48939, 70.83368, 67.882271, 66.51608, 80.4846, 91.017113, 65.768654, 55.468857, 67.30185, 62.439426, 79.17317, 78.26147, 62.73785, 68.438057, 65.845314, 70.724159, 69.812454, 61.86721, 78.02875, 67.72895, 52.57495, 66.08898, 35.501713, 77.41273, 58.75154, 74.45312, 82.302536, 62.398357, 82.91581, 71.71526, 76.22724, 77.17728, 86.15469, 84.665298, 67.742645, 62.73238, 86.264206, 74.5243, 61.289528, 79.983574, 63.50171, 63.501713, 44.38604, 69.37167, 73.508553, 80.958244, 69.64819, 79.786446, 57.30048, 70.313484, 66.8063, 88.982887, 73.724846, 74.81998, 77.697464, 61.163586, 78.962357, 61.02669, 72.72827, 66.004105, 64.11225, 80.216293, 77.719368, 85.23203, 77.577003, 83.230667, 61.77139, 72.32307, 66.03696, 81.497604, 85.801506, 75.61123, 55.468857, 70.9651, 65.700203, 81.916496, 67.589325, 70.1191, 66.13826, 80.95551, 68.3039, 73.163589, 74.685829, 70.65846, 51.55373, 52.6872, 41.5989, 58.13553, 59.16222, 72.780289, 64.991104, 66.75975, 41.20192, 74.27241, 61.93566, 71.504448, 71.485283, 42.962353, 60.64339, 37.52772, 29.70568, 67.644081, 70.31348, 76.648872, 69.93292, 67.04449, 57.9165, 66.6037, 52.67077, 88.23272, 58.49966, 61.275837, 60.7666, 83.720741, 64.142365, 56.654346, 80.106773, 82.51061, 73.127991, 66.86927, 61.68652, 68.77481, 69.577003, 62.62013, 80.810402, 46.97057, 62.0397, 74.74058, 60.64887, 88.164268, 65.75223, 93.149895, 55.93155, 67.78645, 53.864475, 79.304588, 66.036964, 91.84394, 27.23066, 75.42779, 15.3922, 65.39083, 66.77344, 51.16495, 75.28268, 44.79945, 32.7091, 74.05065, 56.19712, 70.18481, 68.29842, 58.15743, 68.928131, 76.02464, 61.0924, 81.57153, 73.60986, 79.90417, 73.17454, 68.8679, 75.70978, 34.77892, 60.91444, 64.449005, 69.05407, 46.96235, 72.01095, 55.30185, 78.63381, 64.06571, 89.05407, 58.42026, 83.98357, 69.85352, 81.34428, 71.35113, 61.85626, 47.36208, 59.37851, 82.46133, 83.73991, 80.17248, 77.74127, 63.15127, 72.64339, 52.61602, 53.53593, 63.50171, 77.226555, 46.63381, 70.598221, 59.115673, 71.88227, 66.47775,
                        60.90075, 75.811089, 67.42231, 65.752228, 77.34154, 58.28337, 55.64955, 67.4798, 69.47296, 48.21629, 60.36687, 76.84326, 80.69267, 61.31417, 61.46475, 78.25599, 67.29911, 54.28337, 31.55921, 65.84257, 65.149895, 67.01163, 53.108829, 68.15058, 70.02053, 57.57974, 67.07461, 69.89733, 64.5202, 65.20192, 49.76591, 79.16222, 74.84463, 73.64271, 52.33128, 65.70842, 72.5065, 67.64134, 79.81109, 73.65366, 68.85421, 65.11978, 68.40794, 80.98289, 84.43532, 88.28474, 79.07735, 68.57495, 74.22314, 93.81519, 52.846, 67.34292, 60.31759, 73.53593, 65.333336, 80.870636, 79.164955, 67.8193, 68.77207, 60.17249, 69.9384, 71.35661, 47.26078, 74.89665, 69.28405, 72.20808, 53.63997, 79.76728, 64.361397, 51.29637, 77.470222, 75.4716, 75.83573, 33.04312, 56.39425, 64.33128, 67.08556, 63.99179, 70.36824, 69.36345, 62.72416, 54.12183, 74.92403, 63.91787, 63.73443, 67.54004, 65.6564, 66.95414, 84.26831, 88.86243, 56.74196, 60.77755, 45.76044, 65.77686, 87.89323, 73.2512, 60.27105, 60.32854, 64.06844, 73.34154, 65.07323, 72.10951, 70.95414, 57.14169, 61.37166, 58.28611, 57.55784, 61.2512, 74.294319, 82.236824, 63.97262, 79.39768, 53.67556, 86.001366, 70.67488, 67.89323, 79.44421, 58.92129, 63.71252, 75.0089, 71.06365, 63.2909, 70.68583, 59.11567, 70.58453, 76.34497, 64.00274, 59.21698, 72.09583, 65.4757, 73.960304, 65.86448, 66.98699, 73.05955, 75.52087, 64.41067, 77.44833, 74.12183, 53.72758, 73.963036, 76.41889, 77.77412, 68.731, 73.96304, 73.01027, 76.547569, 78.58453, 66.19028, 63.17865, 48.88706, 72.43532, 53.80972, 64.72827, 69.00479, 56.08213, 43.718, 62.15195, 65.38809, 62.31622, 69.40178, 75.285423, 78.78987, 63.85216, 68.25736, 86.63381, 84.60232, 86.86379, 70.88844, 61.1499, 69.48939, 63.56194, 58.82546, 58.44216, 66.19849, 44.4271, 39.92334, 73.70294, 70.46133, 78.22588, 56.8898, 72.26283, 65.44285, 58.223133, 72.6078, 53.80698, 74.70225, 80.19165, 28.13142, 67.20055, 78.92677, 62.85558, 64.47639, 66.82272, 71.01163, 77.71937, 44.65709, 67.15127, 63.27173, 72.66804, 62.63929, 60.947296, 54.43943, 71.84942, 64.84873, 80.804932, 63.59753, 66.19028, 50.86105, 81.519508, 72.47091, 78.18481, 60.9473, 80.5859, 61.95483, 65.560577, 69.85078, 62.73511, 63.24983, 74.97878, 84.23546, 71.69336, 82.58453, 64.54483, 61.31417, 81.05407, 62.48049, 70.08624, 78.614647, 71.12389, 67.56742, 77.96304, 75.16222, 65.66735, 65.48666, 87.34018, 65.61259, 31.26078, 63.4935, 76.0, 59.0089, 45.67556, 77.59343, 66.5462, 59.66325, 72.10951, 55.48528, 65.80151, 67.47707, 75.698837, 64.45722, 25.19097, 79.91512, 83.159477, 65.19644, 65.61533, 64.87064, 74.55441, 67.76728, 79.720741, 64.64066, 66.21492, 66.36003, 74.49966, 84.90075, 80.271049, 71.70431, 72.8104, 68.016426,
                        66.27241, 58.35181, 71.83847, 61.66461, 77.05955, 70.62833)
list_of_ages <- list()
starting_index = 59*0 + 1
ending_index = 59*0+ 59 
for (m in starting_index:ending_index) {
  list_of_ages <- list.append(list_of_ages, list_of_ages_all[[m]])
}
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

cat(paste("largest_seen_driver_vaf","        ","largest_neutral_vaf","        ", "cooccurrence"),file="simulated_correlation_plot_Bolton_ages_logconvex_p2_unseen_drivers_first_59.txt",sep="\n")

unseen_fold = 5
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
# u_ben=2*3.8*10^{-5}*3
# u_neu=2*7.4*10^{-5}
u_ben=7.5*10^{-5}*(unseen_fold+1)
u_neu=0.0006
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

          result[[m]]<- mutate(Bio_system[[m]], dt, u, last_mutation_ID,ratio_neutral_beneficial,benefit) #consider the first clone
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
cat(paste(largest_seen_driver_vaf, "        ",largest_neutral_vaf,"        ", cooccurrence_with_seen_drivers),file="simulated_correlation_plot_Bolton_ages_logconvex_p2_unseen_drivers_first_59.txt", append=TRUE, sep="\n")
cooccurrence = 'False'

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







