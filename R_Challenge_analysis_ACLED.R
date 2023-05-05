## R_Challenge_analysis_ACLED.R

"Syd Diller, 4/18/23"

##########################################################################
## INSTRUCTOR: Christopher Fariss
##
## COURSE NAME: Data Science for International Studies (DSIS)
## University of Michigan, Winter 2023, Winter 2022, Winter 2021, Winter 2020
##
## COURSE NAME: Advanced Computational Methods for Social Media and Textual Data (2F)
## University of Essex Summer School 2022, 2021, 2020
##
## Date: 2023-01-20
##
## Please e-mail me if you find any errors or have and suggestions (either email is fine)
## e-mail: cjf0006@gmail.com
## e-mail: cjfariss@umich.edu
##########################################################################
## Instructions:
##
## These challenges are meant to be just that, challenging. They should also be fun. I encourage you to think creatively and collaboratively. Getting stuck or not finishing all the steps is expected and encouraged. This is how learning works.
##
## Always start with step (1) and then continue to each step as time permits.
## Don't worry about completing each step. Document your code for each step.
## You may wish to come back to some of the harder steps as you progress through the course.
## Note that some of the steps may ask you to use skills we have not yet covered in the course.
## Don't worry about these steps now but definitely think through the programming logic if you are stuck and make plans to come back to try them once you feel ready.
##
##########################################################################
##
## Steps for the Challenge
##
## Extract event data from ACLED API at: https://www.acleddata.com/
## Note that because of changes to the ACLED website, I've now provided 3 of the curated datasets from ACLED in the Canvas course website.
## You're welcome to register at ACLED and download additional datasets if you are interested.
##
## (1) Load your ACLED dataset into R.

setwd("C:/Users/sydne/Desktop/R Files/373")

acled_data_LA <- read.csv("LatinAmerica_2018-2023_Feb24.csv")

###########################################################################################################################################
## (2) Inspect and explore the ACLED data.

str(acled_data_LA) ### 255,376 obs. of  29 variables

head(acled_data_LA)

###########################################################################################################################################
## (2a) Specifically, how many rows and columns are contained in the dataset?

nrow(acled_data_LA) ### 255,376

ncol(acled_data_LA) ### 29

###########################################################################################################################################
## (2b) What variable type are each of the variables?

str(acled_data_LA)

"ISO"              ## integer
"EVENT_ID_CNTY"    ## character
"EVENT_ID_NO_CNTY" ## integer
"EVENT_DATE"       ## character
"YEAR"             ## integer
"TIME_PRECISION"   ## integer
"EVENT_TYPE"       ## character
"SUB_EVENT_TYPE"   ## character
"ACTOR1"           ## character
"ASSOC_ACTOR_1"    ## character
"INTER1"           ## integer
"ACTOR2"           ## character
"ASSOC_ACTOR_2"    ## character
"INTER2"           ## integer
"INTERACTION"      ## integer
"REGION"           ## character
"COUNTRY"          ## character
"ADMIN1"           ## character
"ADMIN2"           ## character
"ADMIN3"           ## character
"LOCATION"         ## character
"LATITUDE"         ## numerical
"LONGITUDE"        ## numerical
"GEO_PRECISION"    ## integer
"SOURCE"           ## character
"SOURCE_SCALE"     ## character
"NOTES"            ## character
"FATALITIES"       ## integer
"TIMESTAMP"        ## integer

###########################################################################################################################################
## (2c) What are some the descriptive statistics for the variables (i.e. summary())?

summary(acled_data_LA)

  ### Mean fatalities is only 0.4417, median is 0.0000
  ### GEO_PRECISION median is 1 and mean is 1.141, meaning that the location data is very accurate "If the report notes a particular town, and coordinates are available for that town, the highest precision level “1” is recorded."
  ### TIME_PRECISION median is 1 and mean is 1.133, meaning that most sources include the actual date of the event

###########################################################################################################################################
## (3a) subset() the data to focus on 1 country (use your country from the case selection essay)

Colombia_data <- subset(acled_data_LA, COUNTRY=="Colombia") # subset function to subset acled_data_LA by COUNTRY=="Colombia"
Colombia_data

summary(Colombia_data)
str(Colombia_data)

    ### just for fun / a project for another class

        locations_with_twentyplusfatalities <- Colombia_data$LOCATION[which(Colombia_data$FATALITIES>=10)] ## find locations with 10+ fatalities by accessing corrdinate system 
        locations_with_twentyplusfatalities
        dates_with_twentyplusfatalities <- Colombia_data$EVENT_DATE[which(Colombia_data$FATALITIES>=10)] ## same with dates 
        dates_with_twentyplusfatalities
        
        table(locations_with_tenplusfatalities, dates_with_tenplusfatalities) ### not very useful, just checks that everything lines up 
        paste(locations_with_tenplusfatalities, dates_with_tenplusfatalities, "Fatalities =", Colombia_data$FATALITIES[which(Colombia_data$FATALITIES>=10)], Colombia_data$NOTES[which(Colombia_data$FATALITIES>=10)]) 
        ### prints out the info I need to do more research on specific events 

###########################################################################################################################################
## (3b) Inspect and explore the subset of the original ACLED data

 ### This is code I used to sort out events with the most fatalities and sort it in different ways:

          Colombia_data_highlights <- data.frame(Colombia_data$EVENT_DATE, Colombia_data$YEAR, Colombia_data$EVENT_TYPE, Colombia_data$SUB_EVENT_TYPE, Colombia_data$FATALITIES)
          Colombia_data_highlights ### creates data frame with the variables I am interested in 
          
          sorted_by_FATALITIES <- Colombia_data_highlights[order(Colombia_data_highlights$Colombia_data.FATALITIES, decreasing=TRUE),]  ### Shows data by number of casualties, decreasing
          sorted_by_FATALITIES ### ^^^ Orders the data frame by the number of fatalities of a given event, uses order function within the row section of the dataframe's coordinate system
          
          ### Shows data by year of the event, recent first
          
          ten_or_more_FATALITIES <- which(sorted_by_FATALITIES$Colombia_data.FATALITIES>=10) ## uses which function to select only those entries with more than 10 fatalities
          
          FATALITIES_gtTen_selection <- sorted_by_FATALITIES[ten_or_more_FATALITIES,] ## effectively subsets the events with 10 or more fatalities, just without using the subset fucntion
          
          sorted_by_YEAR <- FATALITIES_gtTen_selection[order(FATALITIES_gtTen_selection$Colombia_data.YEAR, decreasing = TRUE),] ### Shows data by year of the event, recent first
          sorted_by_YEAR ### Year is more useful than Date because the dates aren't in a recognizable format. I could come back transform them into a standard format and then sort, but it is finals week and I am lazy
          ### I did transform dates for a later problem, but don't want to make a whole new data frame for this problem when it is good enough as-is
          
 ### The following is identifying where the most fatalities were in all of Latin America
         allofit <- data.frame(acled_data_LA$LOCATION, acled_data_LA$COUNTRY, acled_data_LA$EVENT_DATE, acled_data_LA$FATALITIES) ### new dataframe with variables that interest me
         Highest_allofit <- allofit[order(allofit$acled_data_LA.FATALITIES, decreasing = TRUE),] ### order by number of fatalities using coordinate system of the df
         Highest_allofit
         
###########################################################################################################################################         
## (4) What is the proportion of EVENT_TYPE for your country?
          
tab1 <- table(Colombia_data$EVENT_TYPE) ### table function shows the frequency of each event type in this data subset
tab1    
# Battles Explosions/Remote   violence                   Protests             Riots   Strategic developments      Violence against civilians 
# 2991                         1062                       5689                 1542         1905                       6048
      

barplot(tab1[order(tab1)], horiz = TRUE, las = 2) ### graphs the frequencies found in tab1. Additional par make the barplot sideways and the [order(tab1)] argument orders the data in tab1 so that it appears in the graph in decending order from most events to least events

proportions(tab1, margin = NULL) ### uses R's proportions function to calculate the proportion of total events that each event type makes up

Violence_against_civillians_Colombia <- table(Colombia_data$EVENT_TYPE == "Violence against civilians") ### code within table fucntion returns logical data for whether an event is violence against civillians (VAC) or not. This is then put into a table to show the total number of VAC vs other event types
proportions(Violence_against_civillians_Colombia, margin=NULL) ### changes the table to show proportions instead of a raw count of event frequency

percent_VaCinColombia <- (((Violence_against_civillians_Colombia[2])*100) / (Violence_against_civillians_Colombia[1]+Violence_against_civillians_Colombia[2])) ### long way of calcualting the % of events that are VAC
aspercent_ColombiaVAC <- paste((percent_VaCinColombia), "% of the events documented in Colombia since 2018 are classifed as Violence Against Civillians") ### paste fucntion used to spit results out in an easy to interpret way
aspercent_ColombiaVAC

#############################################################################################################################################
## (5) Pick one common event in your country. What is the proportion of this EVENT_TYPE in all of the other countries in the dataset?

Violence_against_civillians_regional <- table(acled_data_LA$EVENT_TYPE=='Violence against civilians') ### same code as above but for the whole dataset
proportions(Violence_against_civillians_regional)
    aspercent_regionalVAC <- paste(((((Violence_against_civillians_regional[2]*100) / (Violence_against_civillians_regional[1]+Violence_against_civillians_regional[2])))), "% of the events documented in Colombia since 2018 are classifed as Violence Against Civillians") ### paste function used to spit results out in an easy to interpret way
    aspercent_regionalVAC ###                     ^^^ this could be replaced by proportions(Violence_against_civillians_regional)[2]*100 but I wanted to do it both ways
    
### For each country [individually] in the dataset
    all_countries <- unique(acled_data_LA$COUNTRY) ## use unique function to create vector of every country included in the dataset
    all_countries

  for(i in 1:length(all_countries)){ ### This for loop finds the % of events that are Violence Against Civillians (VAC) for every country in the dataset (when applicable, some countries have no VAC entries)
    
    country_subset <- subset(acled_data_LA, COUNTRY==all_countries[i]) ### use the coordinate syetem of the country name vector all_countries to subset acled_data_LA for each country 
    CountrySub_Table <- table(country_subset$EVENT_TYPE=='Violence against civilians') ### create a table showing the frequency of VAC vs all event types 
    proportion_for_types <- proportions(CountrySub_Table, margin = NULL) ### finds the proportions of VAC vs other types, using the table we just created
    print(paste(all_countries[i], ",", "VAC =", proportion_for_types[2]*100, "%")) ### uses paste function to print the proportion of VAC (converted to %)
### note to self: ***** come back and have this data print into a data frame if you have more time****
  }
    
###########################################################################################################################################
## (6) For just your country, what is the proportion of each of the EVENT_TYPE across the three curated datasets?

gender <- read.csv("gender_Mar19.csv")
    str(gender)
    
    Col_gender <- subset(gender, gender$COUNTRY == 'Colombia') ### subset gender dataset by country of focus
    Col_gen_tab <- table(Col_gender$EVENT_TYPE) # create table of frequencies for event types
    gender_prop <- proportions(Col_gen_tab, margin = NULL) ## convert to proportions
    gender_prop ### use the same code as in prior problem to identify proportions
    
healthworkers <- read.csv("healthworkers_Mar19.csv")
    str(healthworkers)
    
    Col_health <- subset(healthworkers, healthworkers$COUNTRY == 'Colombia')### subset healthworkers dataset by country of focus
    Col_health_tab <- table(Col_health$EVENT_TYPE)
    health_prop <- proportions(Col_health_tab, margin = NULL)
    health_prop
    
journalists <- read.csv("journalists_Mar19.csv")
    str(journalists)
    
    Col_journalists <- subset(journalists, journalists$COUNTRY == 'Colombia')### subset journalists dataset by country of focus
    Col_journalists_tab <- table(Col_journalists$EVENT_TYPE)
    journalists_prop <- proportions(Col_journalists_tab, margin = NULL)
    journalists_prop
    
    #### all of this code is the same as prior problems, just using the new datasets, so see prior problems for more explanation
    
#############################################################################################################################################
## (7a) Calculate the variance of your EVENT_TYPE
    
    ### Formula: [variance] s^2 = sum(term-mean)^2 / (n-1) [n=data points]

  ### long way example: (largely to refresh my own Stats knowledge)
    
            variance_tester <- c(1,7,4,8,3,23,6,4,16,17,19,11,13,19)

            Variance_Function <- function(input){
            
                  themean <- mean(input) ## Step 1: calculate the mean
            
                  result_var_test <- c() 
            
               for(i in 1:length(input)){
                 
                result_var_test[i] <- (input[i] - themean)^2 ## Step 2: subtract the mean and square each integer in the vector
                }
                  print(result_var_test)
                  
                  variance_output <- sum(result_var_test)/(length(input)-1) ## Step 3: sum the numerator and divide by n-1
                  print(paste("variance =", variance_output))
            }
            
            Variance_Function(variance_tester)
            
## finding for my event type
          
  VAC_asBinary <- ifelse(Colombia_data$EVENT_TYPE == "Violence against civilians", 1, 0) ### assign a 1 if an event is VAC and a 0 if not 
  
  var(VAC_asBinary) ### use the binary data we just created to find the variance of VAC vs not VAC

###########################################################################################################################################
## (7b) Calculate the standard deviation of your EVENT_TYPE

  sd(VAC_asBinary) ### use the sd fucntion to find the standard deviation 
  (VAC_asBinary)
            
#############################################################################################################################################
## (8a) Calculate the correlation coefficient for the EVENT_TYPE you selected with one of the other variables or other events in the dataset

        ### Event type by Actor (as integer)
            
  Actorcodes_vs_EventTypes <- table(Colombia_data$INTER1, Colombia_data$EVENT_TYPE) ### make table of actor codes vs event types
        
  for(i in 1:length(unique(Colombia_data$INTER1))){ #### For loop prints the names of all the actors included under each numeric code in Colombia_data$INTER1, use as a key 
    print(i)
    print(unique(Colombia_data$ACTOR1[which(Colombia_data$INTER1==i)]))
      }
  # Inter Code 1: State Forces 
  # Inter Code 2: Rebel Groups 
  # Inter Code 3: Political Militias
  # Inter Code 4: Identity Militias
  # Inter Code 5: Rioters 
  # Inter Code 6: Protesters 
  # Inter Code 7: Civilians 
  # Inter Code 8: External/Other Forces
  
      
EVENTTYPE_numeric <- as.numeric(as.factor(Colombia_data$EVENT_TYPE)) #### changes the names of Event Types into numeric data so that we can calculate variance 
      table(as.numeric(as.factor(Colombia_data$EVENT_TYPE))) ### tables with proportions of event types, using the new numeric codes and the original names
      table(Colombia_data$EVENT_TYPE)
      
      table(EVENTTYPE_numeric)
      
  summary(lm(EVENTTYPE_numeric~Colombia_data$YEAR)) ### combines the summary and lm fucntions to run a linear regression seeing if there is a correlation between the year and the type of events

  cor(x=Colombia_data$YEAR, y=VAC_asBinary) ## use the cor function to find the correlation coefficient of VSC as a binary variable and YEAR 
  cor(x=Colombia_data$YEAR, y=EVENTTYPE_numeric) ### use the cor function to find the correlation coefficient for the numerically coded event types and YEAR                                                          
###########################################################################################################################################
## (8b) Estimate a linear model using your EVENT_TYPE as the left-hand-side variable with one of the other variables or other events in the dataset

  plot(VAC_asBinary~Colombia_data$YEAR, ## Basic plot of Y=VAC_asBinary and X=Colombia_data$YEAR
       data = Colombia_data,
       main = "Event Type over time",
       xlab = "Year",
       ylab = "EVENT TYPE"
      )
  
  regression_data <- lm(VAC_asBinary~Colombia_data$YEAR) ### run linear regression model for Y=VAC_asBinary and X=Colombia_data$YEAR
  summary(regression_data) ### Print Summary 
  
  abline(reg=regression_data) #### plot a line based on the linear regression model 
  ### Not a strong correlation, but seems like VAC has been declining mildly, this lines up with the -0.04919799 correlation coefficient identified above 

###########################################################################################################################################
## (8c) Estimate a linear model using your EVENT_TYPE as the left-hand-side variable with several variables

  str(Colombia_data) ## inspect to select variables 
  
  ### EVENT TYPE and FATALITIES
  
  plot(VAC_asBinary~Colombia_data$FATALITIES, ### basic plot where Y=VAC_asBinary and X=FATALITIES
       data = Colombia_data,
       main = "Event Type vs Fatalities",
       xlab = "Fatalities",
       ylab = "EVENT TYPE"
  )
  
  fit <- lm(VAC_asBinary~Colombia_data$FATALITIES) ### Just like before, run a linear regression and create an object 
  summary(fit) ## r^2 = 0.1515
  cor(x=Colombia_data$FATALITIES, y=VAC_asBinary) ## r = 0.3892266
  
  abline(reg=fit) ### fairly strong positive correlation between VAC event types and the number of fatalities, which makes sense
  
  ### EVENT TYPE and Actor (INTER1)
  
  plot(VAC_asBinary~Colombia_data$INTER1, ## same annotations as above except Y=VAC_asBinary, X=INTER1 (actor codes)
       data = Colombia_data,
       main = "Event Typer vs Actor",
       xlab = "Actor Code",
       ylab = "EVENT TYPE"
  )
  
  fit <- lm(VAC_asBinary~Colombia_data$INTER1)
  summary(fit) ## r^2 = 0.1527
  cor(x=Colombia_data$INTER1, y=VAC_asBinary) # r = -0.390771
  
  abline(reg=fit) ### strong negative correlation between VAC and actor codes, which makes sense
  # Inter Code 1: State Forces 
  # Inter Code 2: Rebel Groups 
  # Inter Code 3: Political Militias
  # Inter Code 4: Identity Militias
  # Inter Code 5: Rioters 
  # Inter Code 6: Protesters 
  # Inter Code 7: Civilians 
  # Inter Code 8: External/Other Forces
  ### More organized actors (1-4) are more likely to inflict VAC than less organized actors/civilians (5-7)
 
  ### EVENT TYPE and LOCATION
  
  Num_Location <- as.numeric(as.factor(Colombia_data$LOCATION)) ## convert Location names into factors and then convert those to numeric data
  
  plot(VAC_asBinary~Num_Location, ### see prior annotations. X=Num_Location (locations as numeric data) Y=VAC_asBinary
       data = Colombia_data,
       main = "Event Type vs Location",
       xlab = "Location",
       ylab = "EVENT TYPE"
  )
  
  fit <- lm(VAC_asBinary~Num_Location)
  summary(fit)## r^2 = 0.007118
  cor(x=Num_Location, y=VAC_asBinary) # r= 0.08436718 pretty much no correlation (as one would expect)
  
  abline(reg=fit) 
  
  ### Very little correlation because Location obviously isn't a scale. If we sorted the location data into categories like "rural, suburban, urban" there might be more meaningful results
  
###########################################################################################################################################
## Warning: The rest of these steps are more challenging:
## 
## GRAPHS:
## (Note we will cover a little bit more about working with text-as-data in 2 weeks)
  
#############################################################################################################################################
## (9) Make a barplot of the proportion of EVENT_TYPE for your country from each of the three datasets included in this challenge
  
  #### Run problem #6 for Set-Up
  
  gender_prop <- proportions(Col_gen_tab, margin = NULL) ### use prop function to find proportions of each event type, this time for the gender curated dataset 
    VACg_prop <- as.vector(gender_prop[5]) ## create object with value of the VAC proportion found in the first step (0.5503876). as.vector is used to remove it from the table and turn it into a vector 
    Gender_vec <- as.vector(gender_prop) ### creates a vector of the proportions found in step 1, using as.vector function to remove the data from the table 
    
  health_prop <- proportions(Col_health_tab, margin = NULL) ## same thing repeated again
    VACh_prop <- as.vector(health_prop[3])
    health_vec <- as.vector(health_prop)
    
  journalists_prop <- proportions(Col_journalists_tab, margin = NULL) ## and again!
    VACj_prop <- as.vector(journalists_prop[4])
    Journalists_vec <- as.vector(journalists_prop)
  
  dev.off() ### clear the plotting space 
  par(mfrow=c(1,1)) ## set parameters (this really isn't necessary here but I want to keep using the par function to make it stick in my brain)
  
      barplot(c(VACg_prop, VACh_prop, VACj_prop), ## create barplot with one bar for each dataset 
                main = "Proportion of Events Classified as Violence Against Civillians", ## add main title 
                names.arg = c("Gender Dataset", "Health Workers Dataset", "Journalists Dataset"), ### add names
                ) ### Barplot kept simple and boring to save time (I spent too much on the next part)
             
  dev.off() ## reset plotting space 
  par(mfrow=c(1,3)) ## create space for 3 plots next to each other 
  
  ### Doing that cool colored barplot thing 
  
      Gender_vec ## Print Gender_vec and gender_prop as reminders about the data we're working with 
      gender_prop
      gender_color <- c("#d00000", "#ffba08", "#3f88c5", "#032b43", "#136f63") ## create color vector the same length as Gender_vec

      barplot(matrix(Gender_vec), col=c(gender_color), xlim=c(0,5), width=3, xlab = "Gender") ### use barplot function to plot Gender_vec as a matrix, using the gender_color vector to assign the correct colors 
      ### add labels and adjust margins and plot width as desired
      
      health_vec ### same steps as above
      health_prop
      health_color <- c("#ffba08", "#3f88c5", "#136f63")
      
      barplot(matrix(health_vec), col=c(health_color), xlim=c(0,5), width=3, xlab = "Health Workers")
      title(main = "Event Types by Dataset") ### I added the main titile here so that it shows up above the center bar 
      
      Journalists_vec ### same again
      journalists_prop
      journalists_color <- c("#ffba08", "#3f88c5", "#032b43", "#136f63")
      
      barplot(matrix(Journalists_vec), col=c(journalists_color), xlim=c(0,5), width=3, xlab = "Journalists")
      
      legend("right", legend=c("Explosions/Remote violence", "Protests", "Riots", "Strategic developments", "Violence against civilians"), fill= c("#d00000", "#ffba08", "#3f88c5", "#032b43", "#136f63"), bty ="o", bg="white") ### create legend on the right side of the plotting space using a vector of colors and respective names for the event types 
                          ### bty ="o", bg="white" makes a white box background for the legend
#############################################################################################################################################
## (10) How many times does the EVENT_TYPE you selected occur in each year of data (2018 - 2020)?
     
      dev.off()
      
      VAC_cada_ano <- table(Colombia_data$EVENT_TYPE=="Violence against civilians", Colombia_data$YEAR)[2,] ### make a table that shows the frequency of violence against civilians as a logical TRUE or FALSE for each year
                                                                                                            ### then use the coordinate system of the table to only access the 2nd row and turn it into an object "VAC_cada_ano"
      barplot(VAC_cada_ano) ### makes a basic barplot with no fancy labels, but provides a basic visualization of changes over the years
      
#############################################################################################################################################
## (11) Create a time series plot for the frequency of each EVENT_TYPE for all countries each week of the year
                        
trialplz <- data.frame( ### creates a dataframe with the following columns:
                        all_types = acled_data_LA$EVENT_TYPE, 
                        protests = replace(acled_data_LA$EVENT_TYPE, acled_data_LA$EVENT_TYPE!="Protests", NA), ### use replace function to replace any EVENT_TYPE in this column that is not "Protests" with an NA
                        VAC = replace(acled_data_LA$EVENT_TYPE, acled_data_LA$EVENT_TYPE!="Violence against civilians", NA), ### same as previous but with "Violence against civilians"
                        SD = replace(acled_data_LA$EVENT_TYPE, acled_data_LA$EVENT_TYPE!="Strategic developments", NA), ### and again
                        Riots = replace(acled_data_LA$EVENT_TYPE, acled_data_LA$EVENT_TYPE!="Riots", NA), ### again 
                        Explosions_RemoteViolence = replace(acled_data_LA$EVENT_TYPE, acled_data_LA$EVENT_TYPE!="Explosions/Remote violence", NA),### otra vez
                        Battles = replace(acled_data_LA$EVENT_TYPE, acled_data_LA$EVENT_TYPE!="Battles", NA), ### Last one 
                        all_dates = as.Date(acled_data_LA$EVENT_DATE, format = "%e-%B-%Y")### Lists dates, but first transforms them into a more compatible format of year-month-day as numbers, not full month names 
)
     
     trialplz[order(trialplz$all_dates, decreasing = FALSE),] #### orders the above dataframe by date, decreasing order 

     freq_dataframe <- data.frame(### Okay so the idea here is to make a data frame that consolidates all of the entries down to a list of each date since the first in the dataset and then report the number of each type of event that took place on that given day
                                  dates = unique(trialplz$all_dates), ### each date in the dataset
                                  protests = data.frame(table(trialplz$all_dates, trialplz$protests))$Freq, ### creates a table that counts the number of "Protests" events for each date in the trialplz data
                                  VAC = data.frame(table(trialplz$all_dates, trialplz$VAC))$Freq,           ### continued^  then it turns that table into a dataframe and uses $Freq to only select the reported frequencies, omitting the dates
                                  SD = data.frame(table(trialplz$all_dates, trialplz$SD))$Freq,             ### all of these columns are doing the same thing with a different event type 
                                  Riots = data.frame(table(trialplz$all_dates, trialplz$Riots))$Freq,
                                  Explosions_RemoteViolence = data.frame(table(trialplz$all_dates, trialplz$Explosions_RemoteViolence))$Freq,
                                  Battles = data.frame(table(trialplz$all_dates, trialplz$Battles))$Freq
                                )

freq_data_ordered <- freq_dataframe[order(freq_dataframe$dates, decreasing = FALSE),] ### Order that dataframe by dates, decreasing

week_start_pos <- seq(from=1, to=nrow(freq_data_ordered), by=7) ### create a sequence of numbers from 1 to the number of rows in the dataframe we just made, but only every 7th number (because weeks have 7 days)
week_starts <- seq(from=min(freq_data_ordered$dates), to=max(freq_data_ordered$dates), by=7) ### use the sequence function again but this time on the dates themselves, to find the start date of each week 

weekly_sum_df <- data.frame(  ### Create an "empty" dataframe
                              Weeks = week_starts, ## one column with all of the week start dates
                              protests = rep(NA, length(week_starts)), ### the rest of the columns just repeat NA for the length of the week_starts vector 
                              VAC = rep(NA, length(week_starts)),
                              SD = rep(NA, length(week_starts)),
                              Riots = rep(NA, length(week_starts)),
                              Explosions_RemoteViolence = rep(NA, length(week_starts)),
                              Battles = rep(NA, length(week_starts))
                            )

    for(i in 1:length(week_start_pos)){ ### create for loop to enter data into the empty dataframe. will repeat a number of times = to the number of week starts in the week_start_pos and week_starts vectors
    
        week <- NULL ## create empty object every time the loop runs
        week <- freq_data_ordered[week_start_pos[i]:(week_start_pos[i]+6),] #loads data from 7 rows of the freq_data_ordered dataframe 
                                                                            # this works by using the for loop repetition and the coordinate system of week_start_pos to create a new selection of rows starting at each integer in week_start_pos
        weekly_sum_df[i,2:7] <- apply(X=week[,2:7], MARGIN = 2, FUN=sum) ### the apply function is used to sum the contents of columns 2:7 by column (not by row), creating one row that is the sums of each column 
                                                                        ### this new sum row is then sent to the correct, corresponding row in the empty  weekly_sum_df dataframe we made previously 
    }

weekly_sum_df ## success! 

dev.off() ### set up plotting space
par(mfrow=c(1,1))

plot.ts(weekly_sum_df[,2:7], plot.type = c("single"), ### plot columns 2:7 (excludes dates) on a single time series 
        xlab= "Weeks (Start Dates)", ## make labels 
        ylab= "Number of Events",
        main = "Number of ACLED Recorded Events by Week (2018-Present)",
        axes = FALSE, ## I do want axes
        col = c("#d00000", "#ffba08", "#3f88c5", "#032b43", "#136f63", "#90be6d") ## Apply colors to each line in the plot, assigned based on order in the dataframe
        )

axis(1, at=1:length(weekly_sum_df$Weeks), labels=weekly_sum_df$Weeks, las=2, cex.axis = 0.45) ## create bottom axis

label_seq <- seq(from=0, to=max(na.exclude(weekly_sum_df$protests)), by=100) ### use sequence function to create a sequence of numerical labels for the y axis 
axis(2, at=label_seq, labels=label_seq, las=2, cex.axis = 0.8) ### plug sequence into the axis funtion to create y axis 
#^^^ this was just excessive and too time consuming

legend("topright", legend=c("Protests", "Violence against civilians", "Strategic developments", "Riots", "Explosions/Remote violence", "Battles"), fill= c("#d00000", "#ffba08", "#3f88c5", "#032b43", "#136f63", "#90be6d"), bty="n") ### creates legend at top with no background

#############################################################################################################################################      
## TEXT ANALYSIS:
## (Note we will cover a little bit more about working with text-as-data in 1 week)
##
## (12) Turn the event description in the ACLED dataset into a Document-by-Term Matrix (DTM). For the DTM, we let i = 1, ..., 
# N index documents and w = 1, ..., W index the unique terms in the collection of documents. For each of the i documents, 
# we determine the frequency of each of the unique $w$ words. Each of the D_iw entries in a DTM represents  
# the number of times the w-th word appears in the i-th document.

Text_data <- acled_data_LA$NOTES ### put the notes into its own object and inspect
head(Text_data)

textEdit <- function(text.vector){ ### use the text edit fucntion as provided, but removing some lines that are
  ### not applicable to this dataset
  TEXT <- text.vector
  TEXT <- gsub("[[:punct:]]", " ", TEXT) # remove all punctuation
  TEXT <- gsub("[^\x20-\x7F\x0D\x0A]", "", TEXT) # remove all non-ascii characters
  TEXT <- gsub("[0-9]", "", TEXT) # remove numbers
  TEXT <- gsub("^\\s+|\\s+$", "", TEXT) # remove extra leading and trailing whitespace
  TEXT <- tolower(TEXT) # turn all letters lower case
  return(TEXT)
} 

setwd("C:/Users/sydne/Desktop/R Files/example data") ## I dont keep everything in the same file so new wd

stopwords <- read.delim("stopwords_twitter.txt", header=FALSE) ### read in stop words 

## for loop set-up:
Event_Types <- as.vector(unique(acled_data_LA$EVENT_TYPE)) ### use unique fucntion to create a vector of all of the event types 

Top_100 <- list(NA,NA,NA,NA,NA,NA) ### Create a new list that is 6 elements long 

for(i in 1:length(Event_Types)){ ### for loop to send data to the list. 
  
  subsetted <- subset(acled_data_LA, acled_data_LA$EVENT_TYPE == Event_Types[i]) ## each iteration of the loop it will subset a different event type 
  
  edittext <- textEdit(subsetted$NOTES) ## use edit text function 
  
  subset_edited_text <- strsplit(edittext, " ") ## split character data into individual words using the strsplit func
  
  for(j in 1:length(subset_edited_text)){ ## for loop to remove filler words 
    subset_edited_text[[j]] <- subset_edited_text[[j]][! subset_edited_text[[j]] %in% as.character(c(stopwords$V1, "no", "size")) ] ## I added no and size because they were very common and not important 
  }
  
  term_freq <- as.data.frame(table(unlist(subset_edited_text))) ### unlist the data, put it in a table to find term frequency, turn that into a data frame
  
  ordered <- term_freq[order(term_freq$Freq, decreasing = TRUE),] ### use coordinate system and order function to sort the term dataframe from highest to lowest frequency 
  
  Top_100[[i]] <- ordered[1:100,] ### take the first 100 rows for each subsetted Event Type notes dataset and insert it into the Top_100 list object we made before this loop
  
}

Top_100 ### now we have a list of 6 dataframes with 100 rows each 

Doc_added <- lapply(1:length(Top_100), function(i){ ### use lapply to create a new list of dataframes using the dataframes we just created but also assigning a document number 
  value <- data.frame(Doc=i, Top_100[[i]])
  return(value)
})

DTM_result <- do.call("rbind", Doc_added) ### use do.call function to rbind the above list of dataframes into one dataframe

### Key for the numbers in Doc for the new dataframe
for(i in 1:length(Event_Types)){ 
  print(paste(i, "=", Event_Types[i]))
}


#############################################################################################################################
## (13) Which words are most commonly associated with the EVENT_TYPE you selected for the entire dataset?

Top_100[[2]] ### prints the second list element of Top_100, which contains data from the Violence Against Civillians subset

dev.off() ## clear plotting space

barplot(height=Top_100[[2]]$Freq[21:2], axes = FALSE, width=1, names.arg=Top_100[[2]]$Var1[21:2], 
        horiz = TRUE, las=2, cex.names = .5, xlab = "Number of Instances",  
        main = "Most Common Terms - Violence against civilians", cex.axis = .75) ### Creates barplot of the top 20 words, using rows 2:21 to avoid row 1 which contains a count of blank spaces 

axis(1, at=seq(from=0, to=max(Top_100[[2]]$Freq[21:2]), by = 1000), 
     labels=seq(from=0, to=max(Top_100[[2]]$Freq[21:2]), by = 1000), 
     las=2, cex.axis = 0.45) ## create bottom axis that counts from 0 to the max frequency of a term, by 1000

#############################################################################################################################
## (14) Create a barplot of the most frequent words by EVENT_TYPE for the entire dataset

Ordered_Top_100_Results <- DTM_result[order(DTM_result$Freq, decreasing = TRUE),] ## sorts the dataframe made at the end of problem 12
unique(Ordered_Top_100_Results$Var1) ## finds all unique terms in DTM_result

dev.off()
par(mfrow=c(2,3)) ### sets plotting space to hold two rows of three plots
### plot the top 20 terms for each event type, using the same barplot arguments as in problem 13 + the colors arguemnt 
barplot(height=Top_100[[1]]$Freq[21:2], width=1, names.arg=Top_100[[1]]$Var1[21:2], horiz = TRUE, las=2, cex.names = .5, xlab = "Number of Instances",  main = "Most Common Terms - Protests", cex.axis = .75, col=c("red"))
barplot(height=Top_100[[2]]$Freq[21:2], width=1, names.arg=Top_100[[2]]$Var1[21:2], horiz = TRUE, las=2, cex.names = .5, xlab = "Number of Instances",  main = "Most Common Terms - Violence against civilians", cex.axis = .75 , col=c("yellow"))
barplot(height=Top_100[[3]]$Freq[21:2], width=1, names.arg=Top_100[[3]]$Var1[21:2], horiz = TRUE, las=2, cex.names = .5, xlab = "Number of Instances",  main = "Most Common Terms - Strategic developments", cex.axis = .75, col=c( "green"))
barplot(height=Top_100[[4]]$Freq[21:2], width=1, names.arg=Top_100[[4]]$Var1[21:2], horiz = TRUE, las=2, cex.names = .5, xlab = "Number of Instances",  main = "Most Common Terms - Riots", cex.axis = .75, col=c("cyan"))
barplot(height=Top_100[[5]]$Freq[21:2], width=1, names.arg=Top_100[[5]]$Var1[21:2], horiz = TRUE, las=2, cex.names = .5, xlab = "Number of Instances",  main = "Most Common Terms - Explosions/Remote violence", cex.axis = .75, col=c("blue"))
barplot(height=Top_100[[6]]$Freq[21:2], width=1, names.arg=Top_100[[6]]$Var1[21:2], horiz = TRUE, las=2, cex.names = .5, xlab = "Number of Instances",  main = "Most Common Terms - Battles", cex.axis = .75, col=c("magenta"))


### create function to Find Term frequency overall using the same process as in problem 12 but without subsetting 
Freq_full_dataset <- function(acled_data_LA){
  
  edittext <- textEdit(acled_data_LA$NOTES) ## edit text using function we made
  subset_edited_text <- strsplit(edittext, " ") ## turn into individual words like before
  
  for(j in 1:length(subset_edited_text)){ ## remove filler words like above
    subset_edited_text[[j]] <- subset_edited_text[[j]][! subset_edited_text[[j]] %in% as.character(c(stopwords$V1, "no", "size")) ] ### same added arguemnts as above 
  }
  
  term_freq <- as.data.frame(table(unlist(subset_edited_text))) ### same process as above
  ordered <- term_freq[order(term_freq$Freq, decreasing = TRUE),]
  return(ordered)
  
}

Term_Freq_Full_Dataset <- Freq_full_dataset(acled_data_LA) ## call function using acled data as the imput and save as new object Term_Freq_Full_Dataset

Termss <- Term_Freq_Full_Dataset[1:50,1] ### isolate the most common terms by accessing the coordinate system of our new dataframe
##  col 1 is terms and we will only be looking at the first 50 because my computer has limitations that I don't know and dont want to reach
empty_df <- data.frame( ### create empty dataframe with columns for the terms vector we just created, the corresponding frequency of those terms for the entire dataset 
  ## and 6 more columns filled with 50 0s using the sequence replicate function. Use 0s instead of NAs because NAs will prevent the barplot function from working later on 
  Term = Termss,
  Freq_Overall = Term_Freq_Full_Dataset[1:50,2],
  Freq.Doc_1 = rep(0, times = length(Termss)),
  Freq.Doc_2 = rep(0, times = length(Termss)),
  Freq.Doc_3 = rep(0, times = length(Termss)),
  Freq.Doc_4 = rep(0, times = length(Termss)),
  Freq.Doc_5 = rep(0, times = length(Termss)),
  Freq.Doc_6 = rep(0, times = length(Termss))
)

list_of_term_info <- NULL ### create new empty object 

for(i in 1:length(Termss)){ ### use for loop to fill the object as follows
  Term_Locations <- which(DTM_result_50$Var1 == Termss[i]) ### identify which rows in DTM_result_50 contain each of the top 50 terms (each for loop focuses on a new term)
  
  Term_Docs <- DTM_result_50$Doc[Term_Locations] ### find which documents contain the top 50 terms using the term locations we found
  Term_Freqs <- DTM_result_50$Freq[Term_Locations] ### do the same but to find the frequencies associated with each term
  
  list_of_term_info[[i]] <- data.frame(Term_Docs, Term_Freqs) ## create new data frame of each of the frequency of each term in each document it is found in
}                                                               ## a data frame is created usning the for loop for each term and then saved as an element of a list in the empty object list_of_term_info

list_of_term_info
names(list_of_term_info) <- Termss ### use the names vector to apply names to each corresponding list element 


for(i in 1:length(list_of_term_info)){ ### this for loop takes the data in the baove list and sorts it into the empty_df we set up above
  
  data_set <- as.data.frame(list_of_term_info[[i]]) ## turn list element into a dataframe and save as 'data'
  
  for(j in 1:nrow(data_set)){ ## use a second nested for loop to apply the term frequencies from list_of_term_info into empty_df
    
    empty_df[i,data_set$Term_Docs[j]+2] <- data_set$Term_Freqs[j] ## use the coordinate system of empty_df to ensure that the frequencies end up in the correct column (this is why there isa +2, to skip the first two columns and get to the frequency columns
  }
}

full_dataset <- empty_df # save the now filled empty_df

dev.off() ## plotting set up
par(mfrow=c(1,1))

barplot(as.matrix(t(full_dataset[2:50,3:8])), names.arg = full_dataset$Term[2:50],  ### use as.matrix and t functions to reformat the now full dataset. Also select [2:50,3:8] to avoid the first row (its just the blank space count) and select only the frequency columns
        beside=FALSE, col=c("#d00000", "#ffba08", "#3f88c5", "#032b43", "#136f63", "#90be6d"),  
        xlab="Frequency", ylab="Terms", main="Common Term Frequency by Event Type", 
        cex.names = .5, cex.lab=1.25, cex.axis = .5, axes = TRUE, las = 2, space = c(0,0))

legend("topright", legend=Event_Types, fill= c("#d00000", "#ffba08", "#3f88c5", "#032b43", "#136f63", "#90be6d"), bg="white") ### create legend on the right side of the plotting space using a vector of colors and respective names for the event types 
### I ran out of time to get the names and axes to fit perfectly
