# Φτιάχνει τις εικόνες (R) και τα HTML (Quarto).
#
#   make          όλα (εικόνες + HTML)
#   make images   μόνο οι εικόνες
#   make html     μόνο τα HTML
#   make clean    σβήνει όλες τις εικόνες και τα HTML
#
# Ξαναφτιάχνεται μόνο ό,τι άλλαξε: αν αλλάξει ένα script, τρέχει μόνο αυτό,
# και μετά ξαναγίνονται τα HTML.

# Εικόνες ανά script
FUNCTION_PNG := values_riemann.png values_mean_all.png \
                values_mean_median_linear.png values_mean_median_convex.png values_mean_median_concave.png \
                values_mad_steps_linear.png values_mad_steps_convex.png values_mad_steps_concave.png
BRIDGE_PNG   := bridge_weights.png
SPREAD_PNG   := mad_steps_full.png mad_steps_interval.png mad_steps_triangle.png \
                var_steps_full.png var_steps_interval.png
MEDIAN_PNG   := normal_mean_median.png normal_mean_median_interval.png

IMAGES := $(FUNCTION_PNG) $(BRIDGE_PNG) $(SPREAD_PNG) $(MEDIAN_PNG)
HTML   := continuous_descriptive_stats.html median.html normal_interval.html

.PHONY: all images html clean

all: html

images: $(IMAGES)

html: $(HTML)

# Κάθε script βγάζει πολλές εικόνες με μία εκτέλεση (grouped target: &:)
$(FUNCTION_PNG) &: function_values.R
	Rscript function_values.R

$(BRIDGE_PNG) &: bridge_visual.R
	Rscript bridge_visual.R

$(SPREAD_PNG) &: mad_visual.R
	Rscript mad_visual.R

normal_mean_median.png: normal_mean_median.R
	Rscript normal_mean_median.R

normal_mean_median_interval.png: normal_mean_median_interval.R
	Rscript normal_mean_median_interval.R

# Κάθε HTML ξαναγίνεται αν αλλάξει το .qmd του ή κάποια εικόνα
%.html: %.qmd $(IMAGES)
	quarto render $<

clean:
	rm -f $(IMAGES) $(HTML)
