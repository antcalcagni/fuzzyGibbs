require(stringr)

exp2pow_all <- function(expr=NULL,print.out=FALSE,print.iter=FALSE){
  K <- sum(str_detect(str_split(expr,pattern = "")[[1]],pattern = "\\^"))
  out <- expr
  for(k in 1:K){
    out <- .exp2pow(out,print.out = print.iter)  
  }
  if(print.out==TRUE){cat("=====================\n");cat("IN: ",expr,"\n\n"); cat("OU: ",out,"\n");cat("=====================\n")}
  return(out)
}

.exp2pow <- function(expr=NULL,print.out=TRUE){
  
  x <- str_split(expr,pattern = "")[[1]]
  i0 <- min(which(str_detect(x,pattern = "\\^")))
  a <- paste(x[1:(i0-1)],collapse = "")
  b <- paste(x[(i0+1):length(x)],collapse = "")
  pow <- paste0("pow( ",.parsing_a(a),",",.parsing_b(b)," )")  #local pow
  
  iid <- str_locate_all(expr,pattern = .str2regex(paste0(.parsing_a(a),"^",.parsing_b(b))))[[1]]
  z <- str_split(expr,pattern = "")[[1]]
  
  out <- NA
  if(iid[1]==1 && iid[2]==length(z)){
    out <- pow
  }else if(iid[1]==1 && iid[2]<length(z)){
    out <- paste(pow,paste(z[(1+iid[2]):length(z)],collapse = ""),collapse = "")  
  }else if(iid[1]>1 && iid[2]==length(z)){
    out <- paste(paste(z[1:(iid[1]-1)],collapse = ""),pow,collapse = "")  
  }else if(iid[1]>1 && iid[2]<length(z)){
    out <- paste(paste(z[1:(iid[1]-1)],collapse = ""),pow,paste(z[(1+iid[2]):length(z)],collapse = ""),collapse = "")
  }
  
  if(print.out==TRUE){cat("=====================\n");cat("IN: ",expr,"\n\n"); cat("OU: ",out,"\n");cat("=====================\n")}
  return(out)
}


.check_exp_conversion <- function(expr_in=NULL,expr_out=NULL,print.out=TRUE){
  out <- -1; xin <- expr_in; xout <- expr_out
  c1 <- str_extract_all(string = xin,pattern = "[a-zA-Z][0-9]+")[[1]] #case: x10, e2 -- name variables containing letters
  c2 <- str_extract_all(string = xin,pattern = "[\\+\\-*\\^\\/][a-zA-Z][\\+\\-*\\^\\/]")[[1]]
    if(length(c2)>0){c2 <- gsub(x = c2,pattern = "[\\-\\+\\*\\^\\/]",replacement = "")} #case: name variables separated by math operations
  c2b <- str_extract_all(string = xin,pattern = "[\\+\\-*\\^\\/][a-zA-Z]")[[1]]
    if(length(c2b)>0){c2b <- gsub(x = c2b,pattern = "[\\-\\+\\*\\^\\/]",replacement = "")} #case: name variables separated by math operations
  c2c <- str_extract_all(string = xin,pattern = "[a-zA-Z][\\+\\-*\\^\\/]")[[1]]
    if(length(c2c)>0){c2c <- gsub(x = c2c,pattern = "[\\-\\+\\*\\^\\/]",replacement = "")} #case: name variables separated by math operations
  c3a <- str_extract_all(string = xin,pattern = "\\([a-zA-Z]\\)")[[1]]
    if(length(c3a)>0){c3a <- gsub(x = c3a,pattern = "\\(|\\)",replacement = "")}
  c3b <- str_extract_all(string = xin,pattern = "\\([a-zA-Z][0-9]+\\)")[[1]] #case: x10, e2 -- name variables containing letters
    if(length(c3b)>0){c3b <- gsub(x = c3b,pattern = "\\(|\\)",replacement = "")}
  
  x <- str_split(xin,pattern = "")[[1]]; if(x[1]=="-"){x[1] <- "-1*"}; xin <- paste(x,collapse = "")
  
  vars <- unique(c(c1,c2,c2b,c2c,c3a,c3b)); vars <- gsub(x = vars,pattern = "+-.",replacement = ""); vars <- vars[sapply(vars,nchar)>0] 
  vals <- runif(length(vars),0.5,2.5)
  for(i in 1:length(vars)){eval(parse(text = paste0(vars[i], " = ", vals[i])))}
  
  pow <- function(x,y){x^y}
  o1 <- eval(parse(text = xin))
  o2 <- eval(parse(text = xout))
  if(print.out){cat(o1," | ",o2,"\n")}
  if(o1==o2){
    out <- 1
  }
  rm(vars,pow)
  return(out)
}


.str2regex <- function(x){
  z <- strsplit(x,split = "")[[1]]
  z <- gsub(x = z,pattern = "\\(",replacement = "\\\\(")
  z <- gsub(x = z,pattern = "\\)",replacement = "\\\\)")
  z <- gsub(x = z,pattern = "\\+",replacement = "\\\\+")
  z <- gsub(x = z,pattern = "\\-",replacement = "\\\\-")
  z <- gsub(x = z,pattern = "\\*",replacement = "\\\\*")
  z <- gsub(x = z,pattern = "\\/",replacement = "\\\\/")
  z <- gsub(x = z,pattern = "\\^",replacement = "\\\\^")
  z <- paste(z,collapse = "")
  return(z)  
}


.parsing_b <- function(chr=NULL){
  
  out <- NULL
  x0 <- str_split(chr,pattern = "")[[1]]
  x00 <- x0[1]; x01 <- x0[2]
  if(length(x0)==1){ #case ^b
    out <- x00
  }else if(str_detect(x00,pattern = "[A-Za-z0-9]") && str_detect(x01,pattern = "[A-Za-z0-9\\.]")){ #check if the ^b is of the tye 2.0 or 20..
    for(i in 2:length(x0)){
      if(!str_detect(x0[i],pattern = "[A-Za-z0-9\\.]")){
         break
      }
    }
    if(i<length(x0)){i <- i-1}
    out <- paste(x0[1:i],collapse = "")
  }else if(str_detect(x00,pattern = "[0-9]")){
    out <- x00
  }else if(str_detect(x01,pattern = "[+*\\-/\\)]")){ #case ^b° or ^b)
    out <- x00
  }else if(str_detect(x00,pattern = "\\(")){ #cases ^(...*), ^(f(..)°)
    z <- x0
    #cbind(z,1:length(z))
    conta_aperte <- 0
    for(i in 2:length(z)){
      if(z[i]=="("){
        conta_aperte <- conta_aperte+1
      }else if(z[i]==")"){
        conta_aperte <- conta_aperte-1
      }
      if(conta_aperte<0){
        break
      }
    }
    out <- paste(z[1:i],collapse="")
  }else if(str_detect(x00,pattern = "[A-Za-z]")){ #case ^f(..)
    for(j in 2:length(x0)){
      if(!str_detect(x0[j],pattern = "[A-Za-z]")){
        break
      }
    }
    z <- x0
    #cbind(z,1:length(z))
    conta_aperte <- 0
    for(i in (j+1):length(z)){
      if(z[i]=="("){
        conta_aperte <- conta_aperte+1
      }else if(z[i]==")"){
        conta_aperte <- conta_aperte-1
      }
      if(conta_aperte<0){
        break
      }
    }
    out <- paste(z[1:i],collapse="")
  } 
  return(out)
}



.parsing_a <- function(chr=NULL){
  
  out <- NULL
  x0 <- str_split(chr,pattern = "")[[1]]
  x00 <- x0[length(x0)-1]; x01 <- x0[length(x0)]
  if(length(x0)==1){ #case a^
    out <- x0
  }else if(str_detect(x01,pattern = "[A-Za-z0-9]") && str_detect(x00,pattern = "[A-Za-z0-9\\.]")){ #check if the a^ is of the tye 2.0 or 20..
    for(i in length(x0):1){ #this also includes cases of type e1, e4, x100
      if(!str_detect(x0[i],pattern = "[A-Za-z0-9\\.]")){
        break
      }
    }
    if(i>1){i <- i+1}
    out <- paste(x0[i:length(x0)],collapse = "")
  }else if(str_detect(x01,pattern = "[0-9]")){
    out <- x01
  }
  else if(str_detect(x00,pattern = "[+*\\-/]") || str_detect(x00,pattern = "\\(")){ #case °a^ or (a^
    out <- x01
  }else if(str_detect(x01,pattern = "\\)")){ #cases (...°)^, (f(..)°)^
    z <- x0
    #cbind(z,1:length(z))
    conta_aperte <- 0
    for(i in (length(z):1)){
      if(z[i]==")"){
        conta_aperte <- conta_aperte+1
      }else if(z[i]=="("){
        conta_aperte <- conta_aperte-1
      }
      if(conta_aperte==0){
        break
      }
    }
    if(i>1 && str_detect(z[i-1],pattern = "[A-Za-z]")){ #case 4: expression in the form f()^
      for(j in (i-1):1){
        if(!str_detect(z[j],pattern = "[A-Za-z]")){
          break
        }
      }
      if(!str_detect(z[j],pattern = "[A-Za-z]")){j <- j+1}
      out <- paste(z[j:length(z)],collapse="")
    }else{
      out <- paste(z[i:length(z)],collapse="")
    }
    
  }
  return(out)
}

R2Cpp_funConverter <- function(funIN=NULL,namefun=NULL,typefun="double",fun_heading=NULL,rcpp_export=TRUE,fileout=NULL,fileout_append=FALSE,filehead=NULL,filetail=NULL,cpp_finalize=FALSE){
  require(stringr)
  
  if(fileout_append==FALSE){
    out <- readLines(con = paste0(getwd(),"/",filehead)) #heading is needed the first time only
  }else{
    out=""
  }
  
  if(fileout_append==FALSE){ #first time the file is printed out
    if(file.exists(paste0(getwd(),"/",fileout))){file.remove(paste0(getwd(),"/",fileout))}
  }
  
  fY_current <- deparse(funIN,width.cutoff = 500) #width.cutoff = 500 should allows not splitting a line over multiple rows
  xargs <- strsplit(strsplit(fY_current[1],split = "function")[[1]][2],split=",")[[1]]
  xargs <- gsub(x = xargs,pattern = "\\(|\\)",replacement = "")
  
  # heading
  hds <- fun_heading
  for(j in 1:length(hds)){
    xargs[j] <- paste(hds[[j]],names(hds)[j])
  }
  fY_current[1] <- paste0(typefun," ",namefun,"(",paste(xargs,collapse = ","),")")
  
  # removing points before variables
  fY_current <- str_replace_all(fY_current,pattern = "\\.([a-zA-Z])",replacement = "\\1")
  
  # add semicolumns, return statement, and proper assignments 
  fY_current <- stringr::str_replace_all(string = fY_current,pattern = "<-",replacement = "=")
  i0 <- max(which(stringr::str_detect(string = fY_current,pattern = "=")))+1
  fY_current[i0] <- paste0("return ",fY_current[i0])
  
  i_last <- length(fY_current) #correct the problem of splitting last line into many sublines
  if((i_last-i0)>1){ #return line split into two or more lines
    x <- paste(fY_current[i0:(i_last-1)],collapse = "")
    x <- unlist(str_replace(string = x,pattern = "return",replacement = ""))
    x <- str_replace_all(string = x,pattern = " ",replacement = "")
    x <- paste0("return   ",x)
    fY_current <- c(fY_current[-c(i0:i_last)],x,"}")
  }
  
  # detecting patterns with sequences: (eg, pars[1:10], pars[(1+ncol(X):ncol(Z))])
  ptn <- "\\[(\\(?\\s*[A-Za-z\\d+*/()\\s]+\\)?\\s*:\\s*\\(?\\s*[A-Za-z\\d+*/()\\s]+\\)?\\s*)\\]"
  iid <- which(stringr::str_detect(string = fY_current,pattern = ptn))
  fY_current[iid] <- tolower(str_replace_all(string = fY_current[iid], "(?i)(ncol|nrow)\\((\\w+)\\)", replacement = "\\2.size()")) #replacement NCOL(X) or NROW(X) with X.size() if any 
  
  iid <- which(str_detect(fY_current,pattern = "ncol|nrow|NCOL|NROW"))
  for(k in iid){
    x0 <- fY_current[k]
    x0 <- str_replace_all(x0,pattern = "ncol|nrow|NCOL|NROW",replacement = "") #in this context vectors can be row-vectors or column-vectors (no diff btw nrow and ncol)
    x0 <- str_replace_all(x0,pattern = "\\(([A-Za-z]+)\\)",replacement = "\\1.size()")
    fY_current[k] <- x0
  }
  
  # detecting patterns with %*% product and substitute with matrixMult1d(..) syntax
  iid <- which(stringr::str_detect(string = fY_current,pattern = "%*%"))
  #ptn <- "(\\b\\w+\\b)\\s*%\\*%\\s*(.+)"
  ptn <- "(\\b\\w+(?:\\[[^\\]]*\\]|\\([^\\)]*\\))?)\\s*%\\*%\\s*(\\b\\w+(?:\\[[^\\]]*\\]|\\([^\\)]*\\))?)"
  fY_current[iid] <- tolower(stringr::str_replace_all(fY_current[iid], ptn, "matrixMult1d(\\1, \\2)"))
  
  # conver ^ to pow(..)
  iid <- which(stringr::str_detect(string = fY_current,pattern = "\\^"))
  for(k in iid){
    x0 <- str_split(fY_current[k],pattern = "=")[[1]] #because we've already substituted the assignment symbol
    if(length(x0)==1){ #length(x0)==1 is the return statement
      x00 <- str_replace(string = x0[1],pattern = "return",replacement = " ")
      x00 <- str_replace_all(string = x00,pattern = " ",replacement = "") #because we've already applied the same function some lines ago
      fY_current[k] <- paste0("return  ",exp2pow_all(expr = x00,FALSE))
    }else{
      #x00 <- trimws(x0[2])
      x00 <- str_replace_all(string = x0[2],pattern = " ",replacement = "")
      fY_current[k] <- paste(x0[1]," = ",exp2pow_all(expr = x00,FALSE))
    } 
  }
  
  # detecting patterns with sequences and substitute with x_seq(..) syntax
  fY_current <- stringr::str_replace_all(fY_current, "(\\b\\w+\\b)\\[(.+?):(.+?)\\]", "x_seq(\\1, \\2, \\3)")
  iid <- which(str_detect(fY_current,pattern = "x_seq"))
  for(k in iid){
    x0 <- str_split(fY_current[k],pattern = "x_seq")[[1]]
    x00 <- str_split(x0[2],pattern = ",")[[1]];
    x00[2] <- str_replace_all(x00[2],pattern = " ",replacement = "")
    if(str_detect(x00[2],pattern = "^\\d+")){ #only number 
      x00[2] <- as.character(as.numeric(x00[2])-1)
    }else if(str_detect(x00[2],pattern = "\\(\\d+$")){ #(number
      x00[2] <-  paste0("(",as.numeric(str_extract(x00[2],pattern = "\\d+"))-1)
    }else if(str_detect(x00[2],pattern = "\\(\\d+[\\+|\\-]")){ #(numer+ or (number-
      x000 <- str_split(trimws(x00[2]),pattern = "\\+")[[1]];
      x000[1] <- paste0("(",as.numeric(str_extract(x00[2],pattern = "\\d+"))-1)
      x00[2] <- paste(x000,collapse = " +")
      }
    x0[2] <- paste0("x_seq",paste(x00,collapse = ",") )
    fY_current[k] <- paste(x0,collapse = " ")
  }
  
  # detecting R's c() operator and convert it properly
  iid <- which(str_detect(fY_current,pattern = "c\\("))
  fY_current[iid] <- str_replace_all(fY_current[iid],pattern = "(\\bc)\\((.+?)\\)",replacement = "\\1({\\2})")
  
  # add last stuffs
  iid <- which(stringr::str_detect(string = fY_current,pattern = "="))
  i0 <- max(iid)+1
  fY_current[iid] <- paste0(fY_current[iid]," ;")
  fY_current[i0] <- paste0(fY_current[i0]," ;")
  fY_current[iid] <- paste0("auto ",fY_current[iid])
  
  # last changes
  fY_current <- tolower(fY_current)
  x <- str_split(str_replace_all(str_replace(fY_current[i0],pattern = "return",replacement = ""),pattern = " ",replacement = ""),pattern = "")[[1]]
  if(x[1]=="-"){
    x[1] <- "-1*"  
    fY_current[i0] <- paste0("return  ",paste(x,collapse = ""))
  }
  
  fY_current <- str_replace_all(fY_current,pattern = "(?<!\\w)gamma\\(",replacement = "tgamma\\(") #resolving name conflicts with the std::tgamma and other similar functions in boost or cmath
  
  
  if(rcpp_export==TRUE){
    out <- c(out," ","// [[Rcpp::export]]",fY_current)  
  }else{
    out <- c(out," ",fY_current)  
  }
  
  if(fileout_append==TRUE){op="a"}else{op="w"}
  conx = file(paste0(getwd(),"/",fileout),open=op)
  writeLines(text = out,conx,sep = "\n",useBytes = TRUE)
  close(conx)
  
  if(cpp_finalize==TRUE){
    tout <- readLines(con = paste0(getwd(),"/",filetail)) #tail
    hout <- readLines(con = paste0(getwd(),"/",fileout))  #current cpp
    out <- c(hout,tout)
    writeLines(con = paste0(getwd(),"/",fileout),text = out,sep = "\n",useBytes = TRUE)
  }
  
  return(invisible())
}

R2Cpp_grad_lnfy_exec <- function(JH=NULL,fileout=NULL,fileout_append=TRUE){
  out <- c("","","// [[Rcpp::export]]")
  
  pars <- paste(paste0("double ",paste0("b",1:JH)),collapse = ",")
  out <- c(out,paste0("std::vector<double> grad_lnfy_exec(double y,vector<double> x,vector<double> z,double lb,double ub, ",pars," ){"))
  
  out <- c(out,
           "",
           "//Note: we cannot use arrays of function pointers here because of the strange behavior of Deriv::Deriv(), which sometimes return vectors instead of scalars.",
           "std::list<std::vector<double>> grad_eval; //temporary container to store gradient elements",
           "std::vector<double> out_grad_eval; //current container for the output",
           "",
           "//because the function based on Deriv::Deriv() returns sometimes a vector instead of a scalar! So safely store the output"
  )
  
  for(j in 1:JH){
    out <- c(out, paste0("vector<double> xout",j," = {grad_lnfy_internal",j,"(y,x,z,lb,ub,",paste(paste0("b",1:JH),collapse = ","),")}; grad_eval.push_back(xout",j,");"))
  }
  
  out <- c(out,"","for(const auto& vec:grad_eval) {out_grad_eval.insert(out_grad_eval.end(), vec.begin(), vec.end());}")
  
  out <- c(out,"","return out_grad_eval;","","}")
  
  if(fileout_append==TRUE){op="a"}else{op="w"}
  conx = file(paste0(getwd(),"/",fileout),open=op)
  writeLines(text = out,conx,sep = "\n",useBytes = TRUE)
  close(conx)
}

R2Cpp_grad_lnfy_exec_loop <- function(JH=NULL,fileout=NULL,fileout_append=TRUE){
  out <- c("","","// [[Rcpp::export]]")
  
  pars <- paste(paste0("double ",paste0("b",1:JH)),collapse = ",")
  out <- c(out,paste0("std::vector<double> grad_lnfy_exec_loop(vector<double>& y,vector<vector<double>>& X,vector<vector<double>>& Z,double lb,double ub, ",pars," ){"))
  
  out <- c(out,
           "",
           "int n = X.size();",
           "int JH = X[0].size() + Z[0].size();",
           "std::vector<std::vector<double>> Out_grad_eval(n,vector<double>(JH));",
           "std::vector<double> out_grad_eval(JH);",
           
           "for (int i = 0; i < n; i++) {",
             "Out_grad_eval[i] = grad_lnfy_exec(y[i], X[i], Z[i],lb,ub,",paste(paste0("b",1:JH),collapse = ","),");",
             "for (int j = 0; j < JH; j++) {",
               "out_grad_eval[j] += Out_grad_eval[i][j];",
             "}",
           "}",
           "",
           "return out_grad_eval;",
           "",
           "}")
  
  if(fileout_append==TRUE){op="a"}else{op="w"}
  conx = file(paste0(getwd(),"/",fileout),open=op)
  writeLines(text = out,conx,sep = "\n",useBytes = TRUE)
  close(conx)
}

R2Cpp_hess_lnfy_exec <- function(K=NULL,JH=NULL,fileout=NULL,fileout_append=TRUE){
  out <- c("","","// [[Rcpp::export]]")
  
  pars <- paste(paste0("double ",paste0("b",1:JH)),collapse = ",")
  out <- c(out,paste0("std::vector<double> hess_lnfy_exec(double y,vector<double> x,vector<double> z,double lb,double ub, ",pars," ){"))
  
  out <- c(out,
           "",
           "//Note: we cannot use arrays of function pointers here because of the strange behavior of Deriv::Deriv(), which sometimes return vectors instead of scalars.",
           "std::list<std::vector<double>> hess_eval; //temporary container to store gradient elements",
           "std::vector<double> out_hess_eval; //current container for the output",
           "",
           "//because the function based on Deriv::Deriv() returns sometimes a vector instead of a scalar! So safely store the output"
  )
  
  for(j in 1:K){
    out <- c(out, paste0("vector<double> xout",j," = {hess_lnfy_internal",j,"(y,x,z,lb,ub, ",paste(paste0("b",1:JH),collapse = ","),")}; hess_eval.push_back(xout",j,");"))
  }
  
  out <- c(out,"","for(const auto& vec:hess_eval) {out_hess_eval.insert(out_hess_eval.end(), vec.begin(), vec.end());}")
  
  out <- c(out,"","return out_hess_eval;","","}")
  
  if(fileout_append==TRUE){op="a"}else{op="w"}
  conx = file(paste0(getwd(),"/",fileout),open=op)
  writeLines(text = out,conx,sep = "\n",useBytes = TRUE)
  close(conx)
}

R2Cpp_hess_lnfy_exec_loop <- function(JH=NULL,fileout=NULL,fileout_append=TRUE){
  out <- c("","","// [[Rcpp::export]]")
  
  pars <- paste(paste0("double ",paste0("b",1:JH)),collapse = ",")
  out <- c(out,paste0("std::vector<double> hess_lnfy_exec_loop(vector<double>& y,vector<vector<double>>& X,vector<vector<double>>& Z,double lb,double ub, ",pars," ){"))
  
  out <- c(out,
           "",
           "int n = X.size();",
           "int JH = X[0].size() + Z[0].size();",
           "",
           "//std::vector<int> l(JH);",
           "//std::iota(l.begin(), l.end(), 1);",
           "//vector<vector<int>> IJ = expand_grid(l,l,true);",
           "//int K = IJ.size();",
           "int K = (JH*(1+JH))/2;",
           "",
           "std::vector<std::vector<double>> Out_hess_eval(n,vector<double>(K));",
           "std::vector<double> out_hess_eval(K);",
           "for (int i = 0; i < n; i++) {",
           paste0("Out_hess_eval[i] = hess_lnfy_exec(y[i], X[i], Z[i],lb,ub,", paste(paste0("b",1:JH),collapse = ","),");"),
             "for (int j = 0; j < K; j++) {",
               "out_hess_eval[j] += Out_hess_eval[i][j];",
             "}",
           "}",
           "",
           "return out_hess_eval;",
           "",
           "}")
  
  if(fileout_append==TRUE){op="a"}else{op="w"}
  conx = file(paste0(getwd(),"/",fileout),open=op)
  writeLines(text = out,conx,sep = "\n",useBytes = TRUE)
  close(conx)
}

