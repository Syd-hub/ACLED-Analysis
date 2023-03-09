setwd()

acled_data_LA <- read.csv("LatinAmerica_2018-2023_Feb24.csv") ### From https://acleddata.com/latin-america-and-the-caribbean/

Colombia_data <- subset(acled_data_LA, COUNTRY = Colombia)
Colombia_data

Colombia_data_highlights <- data.frame(Colombia_data$EVENT_DATE, Colombia_data$YEAR, Colombia_data$EVENT_TYPE, Colombia_data$SUB_EVENT_TYPE, Colombia_data$FATALITIES)
Colombia_data_highlights
colnames(Colombia_data_highlights) <- c("EVENT_DATE", "YEAR", "EVENT_TYPE", "SUB_EVENT_TYPE", "FATALITIES")

sorted_by_FATALITIES <- Colombia_data_highlights[order(Colombia_data_highlights$"FATALITIES", decreasing=TRUE),]  ### Shows data by number of casualties, decreasing
sorted_by_FATALITIES

 ### Shows data by year of the event, recent first

    ten_or_more_FATALITIES <- which(sorted_by_FATALITIES$FATALITIES>=10) ## 233 are 10 or more
    ten_or_more_FATALITIES
    
    FATALITIES_gtTen_selection <- sorted_by_FATALITIES[ten_or_more_FATALITIES,]
    FATALITIES_gtTen_selection
  
  sorted_by_YEAR <- FATALITIES_gtTen_selection[order(FATALITIES_gtTen_selection$"YEAR", decreasing = TRUE),] ### Shows data by year of the event, recent first
  

