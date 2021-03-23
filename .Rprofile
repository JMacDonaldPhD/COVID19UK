



#local(proj_wd <- getwd(), parent.env(environment()))

# Loads package (also changes the parent of the GlobalEnv to allow binding new variables outside Global Environment)
devtools::load_all(".")

# Copies the working directory specific to the user so it can be used
# seamlessly in files which require changing the working directory
local(proj_wd <- getwd(), parent.env(globalenv()))
