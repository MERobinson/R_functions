# function to generate vectorized pptx
library(rvg)
library(officer)
gen_pptx <- function(plot, file, height = 5, width = 5, left = 1, top = 1) {
  read_pptx() %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = dml(ggobj = plot),
            location = ph_location(height = height, width = width,
                                   left = left, top = top),
            bg = "transparent") %>%
    print(target = file)
}

# function to generate combined static/vectorized ppt
gen_statvect_pptx <- function(static, vectorized, file, height = 5, width = 5, left = 1, top = 1) {
  read_pptx() %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = static,
            location = ph_location(height = height, width = width,
                                   left = left, top = top),
            bg = "transparent") %>%
    ph_with(value = dml(ggobj = vectorized),
            location = ph_location(height = height, width = width,
                                   left = left, top = top),
            bg = "transparent") %>%
    print(target = file)
}