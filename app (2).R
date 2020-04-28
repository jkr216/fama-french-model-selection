library(shiny)
library(tidyverse)
library(DT)
library(MASS)
# library(gganimate)
# library(transformr)
# library(png)
# library(gifski)

ui <- fluidPage(
  titlePanel("Interactive Monte-Carlo Option Pricer"), 
  fluidRow(column(12,"This tool compares 3 simple option pricing models and highlights the risk of choosing one particular model.",tags$b("The simulations may take some time."))), fluidRow(column(12, "Please choose the parameters below.")),


  
  sidebarLayout(
    sidebarPanel(
      sliderInput("nb_points", h5("Number of time points"), min = 10, max = 120, value = 50, step = 10),
      sliderInput("nb_sim", h5("Number of Monte-Carlo simulations"), min = 30, max = 2000, value = 300, step = 50),
      sliderInput("nb_sim_shown", h5("Number of simulations shown"), min = 3, max = 24, value = 6),
      numericInput("S", "Stock value (S):", 50, min = 1, max = 100, step = 5),
      numericInput("K", "Strike (K):", 50, min = 1, max = 200, step = 5), 
      numericInput("B", "Barrier (B):", 60, min = 1, max = 200, step = 5), 
      numericInput("r", "Risk-free rate (r):", 0.01, min = 0.00, max = 0.2, step = 0.01),
      numericInput("vol", "Scale (sigma):", 0.2, min = 0.01, max = 1, step = 0.01),
      numericInput("nu", "Nu (Variance Gamma):", 0.2, min = 0.02, max = 1, step = 0.1),
      numericInput("theta", "Theta (Variance Gamma):", 0.2, min = 0.02, max = 1, step = 0.1),
      numericInput("thetah", "Theta (Heston):", 0.04, min = 0.02, max = 1, step = 0.1),
      numericInput("kappa", "Kappa (Heston):", 1.2, min = 0.02, max = 5, step = 0.02),
      numericInput("ksi", "Ksi (Heston):", 0.2, min = 0.02, max = 1, step = 0.02),
      numericInput("rho", "Rho (Heston):", -0.3, min = -0.9, max = 0.9, step = 0.1)
    ),
    mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Models",   fluidRow(column(12,
                                                          "The three models are defined as follows:", 
                                                          withMathJax(),
                                                          helpText(tags$b("1."), 'The Black-Scholes model, such that the underlying has i.i.d. Gaussian log-returns and is expressed under exponential form as 
                                                                   $$S_t=S_0e^{(r-\\sigma^2/2)t+\\sigma B_t}$$ where \\(B_t\\) is a Brownian motion and \\(\\sigma\\) is the volatility of the underlying (of its log-returns). As is customary, \\(r\\) is the risk-free rate.'),
                                                          helpText(tags$b("2."), 'The Variance-Gamma model, for which the underlying is a drifted Brownian motion subject to time-change:
                                                                   $$S_t=S_0e^{(r+w)t+\\theta \\gamma_t^\\nu + \\sigma B_{\\gamma_t^\\nu}},$$
                                                                   where \\(\\gamma_t^\\nu\\) is a gamma subordinator with unit drift and variance \\(\\nu\\). 
                                                                   The martingale condition imposes \\(w=\\text{log}(1-\\theta \\nu-\\sigma^2 \\nu/2)/\\nu\\). The VG process is a Lévy process.'),
                                                         helpText(tags$b("3."), 'The Heston model, in which the underlying and its variance follow the following SDEs:
                                                                  $$\\frac{dS_t}{S_t}=rdt+\\sqrt{v_t}dB_t^1$$ $$dv_t=\\kappa(\\theta-v_t)dt+\\xi \\sqrt{v_t}dB_t^2$$ where the two Brownian motions have \\(\\rho\\) as correlation. \\(\\theta\\) is the long term variance, \\(\\kappa\\) is the speed of mean-reversion for the variance and \\(\\xi\\) determines the variability of the variance.'),
                                                          "The simulation of VG processes follows ", tags$b("Variance Gamma and Monte-Carlo"), "by Michael C. Fu. All others paths are simulated via simple Euler schemes. We impose a floor of 0.001 on the variance in the Heston simulations.  Maturity is set to T=1 for simplicity.", 
                                                          "The formulae for options are the following: 
                                                            $$\\begin{align} \\text{Vanilla Call} & =e^{-rT}\\mathbb{E}[(S_T-K)_+], \\\\ 
                                                         \\text{Asian Call} & =e^{-rT}\\mathbb{E}\\left[\\left(\\frac{1}{T}\\int_0^TS_udu-K\\right)_+\\right], \\\\ 
                                                         \\text{Lookback-Call}&=e^{-rT}\\mathbb{E}\\left[\\underset{u \\le T}{\\max}(S_u)-S_T\\right], \\\\ 
                                                         \\text{Up&In-Call}&=e^{-rT}\\mathbb{E}\\left[(S_T-K)_+1_{\\left\\{\\underset{s \\le T}{\\max}(S_s)>B \\right\\}}\\right],\\\\
                                                         \\text{Down&In-Call}&=e^{-rT}\\mathbb{E}\\left[(S_T-K)_+1_{\\left\\{\\underset{s \\le T}{\\min}(S_s)<B \\right\\}}\\right].
                                                         \\end{align}$$"
                                                          ))),
                    tabPanel("Plots", plotOutput("plot", height = 300), 
                             plotOutput("plot2", height = 300),
                             plotOutput("plot3", height = 300),
                             plotOutput("plot4", height = 300)),
                    tabPanel("Prices", DT::dataTableOutput("prices"))
        )
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  blackscholes <- function(S, X, rf, Tt, sigma) {
    values <- c(2)
    d1 <- (log(S/X)+(rf+sigma^2/2)*Tt)/(sigma*sqrt(Tt))
    d2 <- d1 - sigma * sqrt(Tt)
    values[1] <- S*pnorm(d1) - X*exp(-rf*Tt)*pnorm(d2)
    values[2] <- X*exp(-rf*Tt) * pnorm(-d2) - S*pnorm(-d1)
    values
  }
  
  data <- reactive({ # Generates the data for BS simulations
    S <- input$S
    K <- input$K
    nb_sim <- input$nb_sim
    nb_points <- input$nb_points
    m <- (input$r - input$vol^2 / 2) / nb_points
    vol <- input$vol / sqrt(nb_points)
    sim <- rnorm(n = nb_points * nb_sim, mean = m, sd = vol)
    sim <- rbind(rep(0, nb_sim), matrix(sim, ncol = nb_sim))
    sim <- S * exp(apply(sim, 2, cumsum))
    return(sim)
  })
  
  data_vg <- reactive({ # Generates the data for VG simulations
      S <- input$S
      K <- input$K
      nu <- input$nu
      theta <- input$theta
      nb_sim <- input$nb_sim
      nb_points <- input$nb_points
      vol <- input$vol
      m <- (input$r + log(1-theta*nu-vol^2*nu/2)/nu) / nb_points
      gauss <- rnorm(n = nb_points * nb_sim, mean = 0, sd = 1)
      gam <- rgamma(n = nb_points * nb_sim, shape = 1/nb_points/nu, scale = nu)
      sim <- rbind(rep(0, nb_sim), matrix(m + theta*gam + vol* gauss *sqrt(gam), ncol = nb_sim))
      sim <- S * exp(apply(sim, 2, cumsum))
      return(sim)
  })
  
  data_heston <- reactive({ # Generates the data for VG simulations
      S <- input$S
      K <- input$K
      theta <- input$thetah
      nb_sim <- input$nb_sim
      nb_points <- input$nb_points
      kappa <- input$kappa
      r <- input$r
      ksi <- input$ksi
      rho <- input$rho
      sim <- matrix(0, nrow = (nb_points+1), ncol = nb_sim)
      v <- matrix(0, nrow = (nb_points+1), ncol = nb_sim)
      sim[1,] <- S
      v[1,] <- theta
      Sig <- cbind(c(1,rho),c(rho,1))
      for(i in 1:nb_sim){
          for(j in 2:(nb_points+1)){
              W <- mvrnorm(1, mu = c(0,0), Sigma = Sig)
              v[j,i] <- v[j-1,i] + kappa * (theta - v[j-1,i]) / nb_points + ksi * sqrt(v[j-1,i]/nb_points) * W[1]
              v[j,i] <- max(v[j,i], 0.001) # Hard-coded lower threshold
              sim[j,i] <- sim[j-1,i]*(1+ r/nb_points + sqrt(v[j,i]/ nb_points) * W[2] )
          }
      }
      res <- list(sim = sim, v = v)
      return(res)
  })
  
  prices <- reactive({ # DF with values
    sim <- data()
    pay <- pmax(sim[nrow(sim), ] - input$K,0)   # Empirical Call Price
    e_call <- round(exp(-input$r)*mean(pay), 3)
    bs_call <- round(blackscholes(input$S, input$K, input$r, 1, input$vol)[1],3)
    if(input$B<input$S & input$B<input$K){
    bs_dic <- round(blackscholes(input$B^2/input$S, input$K, input$r, 1, input$vol)[1]*(input$B/input$S)^(2*input$r/input$vol^2-1),3)}
    else {bs_dic <- bs_call}
    avg <- pmax(apply(sim,2,mean) - input$K,0)  # Asian 
    a_call <- round(exp(-input$r)*mean(avg), 3)
    mx <- apply(sim,2,max)                      # Lookback & Barrier
    cond <- (mx > input$B)
    ui_call <- round(exp(-input$r)*mean(cond*pay),3)
    mi <- apply(sim,2,min)                      # Barrier
    cond <- (mi < input$B)
    di_call <- round(exp(-input$r)*mean(cond*pay),3)
    l_call <- round(exp(-input$r)*mean(mx-input$S), 3)
    
    sim <- data_vg()
    pay <- pmax(sim[nrow(sim), ] - input$K,0)   # Empirical Call Price
    e_call_vg <- round(exp(-input$r)*mean(pay), 3)
    avg <- pmax(apply(sim,2,mean) - input$K,0)  # Asian 
    a_call_vg <- round(exp(-input$r)*mean(avg), 3)
    mx <- apply(sim,2,max)                      # Lookback & Barrier
    cond <- (mx > input$B)
    ui_call_vg <- round(exp(-input$r)*mean(cond*pay),3)
    mi <- apply(sim,2,min)                      # Barrier
    cond <- (mi < input$B)
    di_call_vg <- round(exp(-input$r)*mean(cond*pay),3)
    l_call_vg <- round(exp(-input$r)*mean(mx- input$S), 3)
    
    sim <- data_heston()$sim
    pay <- pmax(sim[nrow(sim), ] - input$K,0)   # Empirical Call Price
    e_call_h <- round(exp(-input$r)*mean(pay), 3)
    avg <- pmax(apply(sim,2,mean) - input$K,0)  # Asian 
    a_call_h <- round(exp(-input$r)*mean(avg), 3)
    mx <- apply(sim,2,max)                      # Lookback & Barrier
    cond <- (mx > input$B)
    ui_call_h <- round(exp(-input$r)*mean(cond*pay),3)
    mi <- apply(sim,2,min)                      # Barrier
    cond <- (mi < input$B)
    di_call_h <- round(exp(-input$r)*mean(cond*pay),3)
    l_call_h <- round(exp(-input$r)*mean(mx- input$S), 3)    
    Prices <- c("Black-Scholes Vanilla Call", "Simulated Vanilla Call", "Simulated Asian Call", 
                "Simulated Lookback Call", "Simulated Up&In Call Price", "Simulated Down&In Call",
                "Black-Scholes Down&In Call")
    BS <- c(bs_call, e_call, a_call, l_call, ui_call, di_call,bs_dic)
    VG <- c(NA, e_call_vg, a_call_vg, l_call_vg, ui_call_vg, di_call_vg,NA)
    Heston <-  c(NA, e_call_h, a_call_h, l_call_h, ui_call_h, di_call_h,NA)
    tab <- data.frame(Prices,BS,VG, Heston)
    colnames(tab) <- c("Prices", "Black-Scholes", "Variance Gamma", "Heston")
    return(tab)
  }) 
  
  output$prices <- DT::renderDataTable({prices()})
  
  output$plot <- renderPlot({ # The plot
    nb_sim <- input$nb_sim
    nb_sim_shown <- input$nb_sim_shown
    nb_points <- input$nb_points
    sim <- data.frame(cbind(0:nb_points, data()))
    sim <- sim[,1:(nb_sim_shown+1)]
    colnames(sim) <- c("Iteration", 1:nb_sim_shown)
    sim <- gather(sim, key = Simulation, value = Value, -Iteration)
    sim$Simulation <- as.factor(sim$Simulation)
    p <- ggplot(sim, aes(x = Iteration, y = Value, color = Simulation)) + geom_line() + theme(text = element_text(size=16)) +
      geom_hline(yintercept = input$K) +  geom_hline(yintercept = input$B, linetype = 2) + 
        annotate("text", 1, input$B, vjust = -1, label = "Barrier") +
      annotate("text", 1, input$K, vjust = -1, label = "Strike") + # + transition_time(as.numeric(Simulation))
      labs(title = "Gaussian simulations")
     print(p)
  }, height = 300)
  
  output$plot2 <- renderPlot({ # The plot
      nb_sim <- input$nb_sim
      nb_sim_shown <- input$nb_sim_shown
      nb_points <- input$nb_points
      sim <- data.frame(cbind(0:nb_points, data_vg()))
      sim <- sim[,1:(nb_sim_shown+1)]
      colnames(sim) <- c("Iteration", 1:nb_sim_shown)
      sim <- gather(sim, key = Simulation, value = Value, -Iteration)
      sim$Simulation <- as.factor(sim$Simulation)
      p <- ggplot(sim, aes(x = Iteration, y = Value, color = Simulation)) + geom_line() + theme(text = element_text(size=16)) +
          geom_hline(yintercept = input$K) + geom_hline(yintercept = input$B, linetype = 2) + 
          annotate("text", 1, input$B, vjust = -1, label = "Barrier") +
          annotate("text", 1, input$K, vjust = -1, label = "Strike") + # + transition_time(Simulation)
          labs(title = "Variance Gamma simulations")
      print(p)
  }, height = 300)
  
  output$plot3 <- renderPlot({ # The plot
      nb_sim <- input$nb_sim
      nb_sim_shown <- input$nb_sim_shown
      nb_points <- input$nb_points
      sim <- data.frame(cbind(0:nb_points, data_heston()$sim))
      sim <- sim[,1:(nb_sim_shown+1)]
      colnames(sim) <- c("Iteration", 1:nb_sim_shown)
      sim <- gather(sim, key = Simulation, value = Value, -Iteration)
      sim$Simulation <- as.factor(sim$Simulation)
      p <- ggplot(sim, aes(x = Iteration, y = Value, color = Simulation)) + geom_line() + theme(text = element_text(size=16)) +
          geom_hline(yintercept = input$K) + geom_hline(yintercept = input$B, linetype = 2) + 
          annotate("text", 1, input$B, vjust = -1, label = "Barrier") +
          annotate("text", 1, input$K, vjust = -1, label = "Strike") + # + transition_time(Simulation)
          labs(title = "Heston simulations")
      print(p)
  }, height = 300)
  
  output$plot4 <- renderPlot({ # The plot
      nb_sim <- input$nb_sim
      nb_sim_shown <- input$nb_sim_shown
      nb_points <- input$nb_points
      sim <- data.frame(cbind(0:nb_points, sqrt(data_heston()$v)))
      sim <- sim[,1:(nb_sim_shown+1)]
      colnames(sim) <- c("Iteration", 1:nb_sim_shown)
      sim <- gather(sim, key = Simulation, value = Value, -Iteration)
      sim$Simulation <- as.factor(sim$Simulation)
      p <- ggplot(sim, aes(x = Iteration, y = Value, color = Simulation)) + geom_line() + theme(text = element_text(size=16)) +
          labs(title = "Heston volatility")
      print(p)
  }, height = 300)
}

# Run the app ----
shinyApp(ui = ui, server = server)
