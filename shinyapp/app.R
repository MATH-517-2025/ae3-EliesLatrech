#Assignment 3 - Eliès Latrech
#Can use this shiny app to observe the effects of the parameters on the estimate of the optimal global bandwidth

library(shiny)
library(ggplot2)

# true function m
m_fun <- function(x) sin(1 / (x/3 + 0.1))

# 2nd derivative from quartic polynomial coefficients
m2_from_coefs <- function(coefs, x) {
  2*coefs[3] + 6*coefs[4]*x + 12*coefs[5]*x^2
}

# estimate sigma2, theta22, h_AMISE for given N
estimate_hamise <- function(X, Y, N) {
  n <- length(X)
  breaks <- seq(min(X), max(X), length.out = N+1)
  block_idx <- cut(X, breaks, include.lowest = TRUE, labels = FALSE)
  
  RSS <- 0
  second_sq <- rep(NA, n)
  
  for (j in 1:N) {
    sel <- which(block_idx == j)
    if (length(sel) >= 5) {
      df <- data.frame(x = X[sel], y = Y[sel])
      fit <- lm(y ~ poly(x, 4, raw = TRUE), data = df)
      coefs <- coef(fit)
      coefs <- c(coefs[1], coefs[2:5] %||% 0)
      
      fitted <- predict(fit, newdata = df)
      RSS <- RSS + sum((df$y - fitted)^2)
      second_sq[sel] <- m2_from_coefs(coefs, df$x)^2
    }
  }
  
  sigma2_hat <- RSS / (n - 5*N)
  theta22_hat <- mean(second_sq, na.rm = TRUE)
  supp_len <- max(X) - min(X)
  
  h_hat <- n^(-1/5) * ((35 * sigma2_hat * supp_len) / theta22_hat)^(1/5)
  
  list(h = h_hat, sigma2 = sigma2_hat, theta22 = theta22_hat, RSS = RSS)
}

# Mallow's statistic
Cp_value <- function(RSS, n, N, RSS_max, Nmax) {
  RSS / (RSS_max / (n - 5*Nmax)) - (n - 10*N)
}

ui <- fluidPage(
  titlePanel("h_AMISE Explorer"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n", "Sample size n:", 50, 1000, 200, step = 10),
      sliderInput("alpha", "Beta alpha:", 0.5, 5, 2, step = 0.1),
      sliderInput("beta", "Beta beta:", 0.5, 5, 2, step = 0.1),
      sliderInput("N", "Blocks N:", 1, 30, 5, step = 1),
      sliderInput("sigma2", "Variance of normal noise:", 0, 1.5, 1, step=0.1),
      actionButton("optN", "Set N to Cp-optimal")
    ),
    mainPanel(
      verbatimTextOutput("vals"),
      plotOutput("plot")
    )
  )
)

server <- function(input, output, session) {
  data_gen <- reactive({
    X <- rbeta(input$n, input$alpha, input$beta)
    Y <- m_fun(X) + rnorm(input$n, sd = input$sigma2)
    list(X = X, Y = Y)
  })
  
  current_N <- reactiveVal(5)
  observe({ current_N(input$N) })
  
  # button to set N optimally
  observeEvent(input$optN, {
    dat <- data_gen()
    n <- input$n
    Nmax <- max(min(floor(n/20), 5), 1)
    RSSmax <- estimate_hamise(dat$X, dat$Y, Nmax)$RSS
    gridN <- 1:Nmax
    Cp_vals <- sapply(gridN, function(N) {
      RSS <- estimate_hamise(dat$X, dat$Y, N)$RSS
      Cp_value(RSS, n, N, RSSmax, Nmax)
    })
    bestN <- gridN[which.min(Cp_vals)]
    updateSliderInput(session, "N", value = bestN)
    current_N(bestN)
  })
  
  output$vals <- renderPrint({
    dat <- data_gen()
    res <- estimate_hamise(dat$X, dat$Y, current_N())
    cat("n =", input$n, " alpha =", input$alpha, " beta =", input$beta,
        " N =", current_N(), " sigma squared =", input$sigma2, "\n\n")
    cat("Estimated sigma^2:", round(res$sigma2, 4), "\n")
    cat("Estimated theta22:", round(res$theta22, 4), "\n")
    cat("h_AMISE:", round(res$h, 4), "\n")
  })
  
  output$plot <- renderPlot({
    dat <- data_gen()
    X <- dat$X; Y <- dat$Y
    N <- current_N()
    breaks <- seq(min(X), max(X), length.out = N+1)
    block_idx <- cut(X, breaks, include.lowest = TRUE, labels = FALSE)
    df <- data.frame(x = X, y = Y, block = block_idx)
    
    ggplot(df, aes(x, y, color = factor(block))) +
      geom_point() +
      labs(title = paste("n =", input$n, ", N =", N),
           color = "Block") +
      theme_minimal()
  })
}

shinyApp(ui, server)
