setwd("C:/Users/sydne/Documents/Winter '23/373")

acled_data_LA <- read.csv("LatinAmerica_2018-2023_Feb24.csv")

Colombia_data <- subset(acled_data_LA, COUNTRY = Colombia)
Colombia_data

Colombia_data_highlights <- data.frame(Colombia_data$EVENT_DATE, Colombia_data$EVENT_TYPE, Colombia_data$SUB_EVENT_TYPE, Colombia_data$FATALITIES)
colnames(Colombia_data_highlights) <- c("EVENT_DATE", "EVENT_TYPE", "SUB_EVENT_TYPE", "FATALITIES")
Colombia_data_highlights

Colombia_data_highlights[order(Colombia_data_highlights$"FATALITIES", decreasing=TRUE),] 
