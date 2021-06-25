#prep data for OSTATs----
library(dplyr)
library(purrr)
library(Ostats)

bird_data<-read.csv("~/Documents/temp/ben_birds_forQ.csv")


#these are the locations that were in the plot. specifically VINS Iis the one that is curious
sub_station<-c("VINS","PATT","FTGI")


dat_in <- bird_data %>%
  filter(STATION %in% sub_station)%>% #if you want to just run the three staions above, if you want to run them all comment out this line (takes ~30 min to run them all)
  select(STATION, SPEC, WEIGHT) %>%
  filter(!is.na(WEIGHT)) %>%
  mutate(log_WEIGHT = log10(WEIGHT))


# Load the modified community_overlap function that returns raw values too
source('community_overlap_modified.R')

# Now do the community overlap for each site, using the modified function that also returns the raw list.
raw_overlap_list <- dat_in %>%
  group_by(STATION) %>%
  group_map(~ community_overlap(traits = as.matrix(.[, 'log_WEIGHT']), sp = factor(.$SPEC), normal = TRUE, output = 'median'))

# We have a list of 3 lists (one for each site). Each site has an element "raw" which is a dataframe with the raw overlaps.
# Combine them into a single dataframe and write to csv.

raw_overlap_df <- map2_dfr(sub_station, raw_overlap_list, ~ data.frame(STATION = .x, .y[['raw']]))

write.csv(raw_overlap_df, 'bird3stationoverlaps.csv', row.names = FALSE)

# Let's take a look at the distribution of overlaps by site

library(ggplot2)

ggplot(raw_overlap_df, aes(x = overlap)) +
  geom_histogram() +
  facet_wrap(~ STATION) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)))

# Proportion very close to zero overlaps, by site
raw_overlap_df %>%
  group_by(STATION) %>%
  summarize(n_verysmall = sum(overlap < 1e-5),
            n_notsmall = sum(overlap > 1e-5),
            proportion_verysmall = n_verysmall/(n_verysmall+n_notsmall))


# try different kernel bandwidth ------------------------------------------

# The larger the bandwidth, the "wider" the distributions are around the data points
# By default, the Ostats function uses the default behavior of R's density() function to choose a bandwidth in every case
# This uses an algorithm called "nrd0".
# We can set a fixed bandwidth value instead.

# Let's try a few.



# Example of one distribution.
chickadee <- filter(dat_in, SPEC == 'CACH', STATION == 'FTGI')

chick_dens <- map(list('nrd0', 0.01, 0.05, 0.10), ~ density(chickadee[, 'log_WEIGHT'], bw = .))

# Plot the densities
chick_dens_df <- map_dfr(chick_dens, ~ data.frame(bw = .[['bw']], x = .[['x']], y = .[['y']]))
ggplot(chick_dens_df, aes(x = x, y = y, color = factor(bw), group = factor(bw))) +
  geom_line(size = 1) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)))

# Now, let's try to do the entire thing with different bandwidths
raw_overlap_list_bandwidths <- map(list('nrd0', 0.01, 0.05, 0.10), function(bandwidth) dat_in %>%
                                     group_by(STATION) %>%
                                     group_map(~ community_overlap(traits = as.matrix(.[, 'log_WEIGHT']), sp = factor(.$SPEC), normal = TRUE, output = 'median', density_args = list(bw = bandwidth)))
)

# Extract the values and look at them in a table.
ostats_list <- map(raw_overlap_list_bandwidths, ~ do.call(rbind, map(., "value")))
ostats_table <- data.frame(Station = sub_station, ostats_list)
names(ostats_table) <- c('Station', 'bw = default', 'bw = 0.01', 'bw = 0.05', 'bw = 0.10')

# Replace very small values with a string and make a table
ostats_table %>%
  mutate(across(where(is.numeric), function(x) ifelse(x < 1e-6, "< 1e-6", round(x, 3)))) %>%
  knitr::kable()

# Looks like default is between 0.01 and 0.05.
