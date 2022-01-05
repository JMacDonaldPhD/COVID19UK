# Disease Progression Graph


# Directed Graph

# 3 status of the connection between two states
# 0 - not connected, i.e. you cannot reach one state to the other (directly)
# 1 - connected, progression is independent of other agents
# 2 - connected, progression is dependent on other agents


# As well as the definition of connections, there are then the equations
# which govern the rate of transition from state to state.

# Each edge must be associated with an rate equation
# The nature of this rate equation depends on the connection between two states
# is of type 1 or type 2 (i.e independent/dependent)

# Rate equations for type 1 connections should be dependent only on an
# individuals characteristics. These chararcteristics may change over time,
# but for simplicity they are assumed to stay constant for now.
# This means the rate of transition between two states of connection type 1
# can be precalculated. If the characteristics can change throughout time,
# then precalculation may or may not waste time over the course of a simulation

# Rate equations for type 2 connections will most likely have to be recalculated
# after an event occurs.
# Detect when a variable which this connection is reliant on changes to decide
# whether calculations have to take place.

# SIR Example


# The below code is a bit annoying to write for a user.

# They would probably like to say there are states S, I and R
# Then state that S goes I and I goes R
# Can the rate equations then infer whether the connections are of type 1 or
# 2 or is it better to leave that for the user to define?


# Idea
x <- c(S + I ~ I, I ~ R)
diseaseGraph <- matrix(c(0, 2, 0,
                         0, 0, 1,
                         0, 0, 0),
                       ncol = 3, nrow = 3, byrow = T)

rownames(diseaseGraph) <- colnames(diseaseGraph) <- c("S", "I", "R")

rateEquations <- list(infection = function(X, beta) X[2]*beta/N,
                      removal = function(gamma) gamma)


infectionRate <- expression(I*beta/N)
# extract Infection Rate dependencies

as.character(infectionRate)
expr_dep <- function(expr){

}
substr("hello", 1, 1)

exp <- grep("A{3}", c("AAA", "AA"), value = T)

# Look for brackets so we don't look in there for
# subtraction




regexpr("\\*", c("*", infectionRate))


brackets_regex <- "\\([^\\)]*\\)"
subtract_regex <- "(\\(.*?\\)[^-]+|[^-]+)"

gregexpr(subtract_regex, string, perl = T)

# Look for subtraction outside brackets
match_subtract <- function(string){

  subtract_regex <- "(\\(.*?\\)[^-]+|[^-]+)"
  sub_match <- gregexpr(subtract_regex, string, perl = T)


  match_indices <- sub_match[[1]][c(1,2, 3)]
  match_lengths <- attributes(sub_match[[1]])$match.length
  print(match_indices)

  left <- substr(string, match_indices[1],
                 match_indices[1] + match_lengths[1] - 1)
  right <- substr(string, match_indices[2],
                  match_indices[2] + match_lengths[2] - 1)


  return(list(left = left, right = right))

}
string <- "(A*B - B) + ((B*A) - A))"

debug(match_subtract)
match_subtract(string)

match_addition <- function(string){
  subtract_regex <- "(\\(.*?\\)[^+]+|[^+]+)"
  sub_match <- gregexpr(subtract_regex, string, perl = T)


  match_indices <- sub_match[[1]][c(1,2, 3)]
  match_lengths <- attributes(sub_match[[1]])$match.length
  print(match_indices)

  left <- substr(string, match_indices[1],
                 match_indices[1] + match_lengths[1] - 1)
  right <- substr(string, match_indices[2],
                  match_indices[2] + match_lengths[2] - 1)


  return(list(left = left, right = right))
}

match_addition(string)
match_subtract(string)
sub_match[[1]][c(1,2)]

"\\+"

"\\*"

"/"

"\\^"


"(" # then find last match of ")"

"[^-]*\([^\)]+\)([^-]+|.*$"



# Use this function to find the first instance of
# a symbol that isn't nested within round brackets
find_unnested_symbol <- function(symbol, string){

  if(nchar(symbol) > 1){
    stop("Provide symbol of length 1")
  }

  nestLevel <- 0
  l <- nchar(string)
  found <- FALSE

  # Find all instances of open and closed brackets
  # and symbol.
  patt <- paste0(c("[\\(\\)", symbol, "]{1}"),
                 collapse = "")
  matches <- gregexpr(patt, string, perl = T)


  noMatches <- length(matches[[1]])
  matchIndices <- matches[[1]][1:noMatches]
  i <- 1

  # If nestLevel == 0, find next open bracket or
  # symbol. If open bracket comes first increase
  # nestLevel by 1. While nestLevel > 0 search for
  # open or closed brackets. If closed comes first
  # decrease nestLevel by 1. If open comes first
  # Increase nest level by one. Continue until
  # symbol is found at nestLevel zero
  while(i <= noMatches & !found){
    match <- substr(string, matchIndices[i],
                    matchIndices[i])
    if(nestLevel == 0){
      if(match == "("){
        nestLevel <- nestLevel + 1
      } else if(match == symbol){
        found = TRUE
      } else if(match == ")"){
        stop("Invalid Equation given")
      }
    } else{
      if(match == "("){
        nestLevel <- nestLevel + 1
      } else if(match == ")"){
        nestLevel <- nestLevel - 1
      }
    }
    i <- i + 1
    #print(nestLevel)
  }
  return(list(found = found, index = matchIndices[i - 1]))
}


# Attempts to find the arguments the provided binary
# operator is acting on (not nested)
binop_arguments <- function(string, op){

  if(!(op %in% c("+", "-", "*", "/"))){
    stop("Invalid Binary Operator provided")
  }

  sub_index <- find_unnested_symbol(op, string)

  if(sub_index$found){
    left <- substr(string, start = 1, stop = sub_index$index - 1)
    right <- substr(string, start = sub_index$index + 1, stop = nchar(string))
    return(list(left = left, right = right))
  } else{
    message("Binary not found operation found. Returning string")
    return(list(left = string))
  }
}

binop_arguments("(A - B) + C", "/")


expression_parse <- function(expr){
  # First check whether expression is in brackets

  # Then check for operators in this order
  # 1. Subtraction
  # 2. Addition
  # 3. Multiplication
  # 4. Division
  # 5. Indicies
  # 6. Brackets
  ops <- c("-", "+", "*", "/", "^", "()")

  for(i in 1:(length(ops) - 1)){
    expr <- binop_arguments(expr, ops[i])

    if(!(is.null(expr$right))){
      break
    } else{
      expr <- expr$left
    }
  }
  return(list(arg = expr, op = ops[i]))
}

address(N)
?address

find_unnested_symbol("-", "(A-B) * A")

expression_parse("(A - B) * A")


# Precedence of operators
# Brackets
# Indicies
# Division
# Multiplication
# Addition
# Subtraction

# More complicated with other operations such as functions

# Assume there are no function, just primative operators
# Look in reverse order (look for things not inside brackets at first)
# So look for subtraction and seperate the two components
# If no substract then do addition




attributes(infectionRate)

X <- c(90, 10, 0)
I <- X[2]
beta <- 1
N = sum(X)
eval(infection)
eval(expression(X[2]*beta/N))



# List rate equation dependencies

args(rateEquations$infection)

attr(rateEquations$infection)

# Sort individuals by characteristics to reduce the number of type of events
# to choose between. I.e there is no need to define which parts of the
# population are homogeneous etc. The model will figure that out based on
# the characteristics provided.

# Is this ideal, or would users like to be specific about whether a population
# is homogeneous?

# Final(?) element is that of mixing assumptions.
# Is there metapopulations?
# Is there household/school/community structures?




postfix_expr <- function(expr){
  char <- as.character(expr)

}


stack <- function(){
  new
}









