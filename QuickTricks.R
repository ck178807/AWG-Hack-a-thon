#example 1
#Function can be used to break up cells in an r data frame by "split".  
#In this case, this would take the 1st field in each cell in the specified column of df and break it apart with comma as the delimiter and return the first field.
sapply(strsplit(df$column,split = ","),'[', 1)

#example 2
#paste function can be used to generate combined lists
#(assuming the metadata is a set of columns that precede or follow the counts/intensities data...See example 4)
#For example, the command below would allow you to add a new column to a dataframe that combined the gravity metadata column with the genotype column
df$Grav_Geno<-paste(df$Gravity, df$Genotype,sep="_")

#example 3
#thinning data frames
#(assuming the metadata is a set of rows stacked with the counts/intensities data...see example 4)
#you could do the following to thin the data frame to only contain the Microarray experiments
df_MAonly<-df[grep(df$Assay,pattern = "Microarray"),]
#and this can be further extended to include only root tissue samples, assayed by microarray whose age is 8 days or less as follows
df_MAonly<-df[grep(df$Assay,pattern = "Microarray") & grep(df$Assay,pattern = "Roots") & df$Age <= 8,]