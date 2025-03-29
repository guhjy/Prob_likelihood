library(shiny)
library(HDInterval)

# --- UI Definition ---
ui <- fluidPage(
  # Update the title to reflect the new explanation framework
  titlePanel("機率與似然率: 頻率 vs 貝氏學派"),

  # Enable MathJax for rendering mathematical formulas
  withMathJax(),

  sidebarLayout(
    # --- Sidebar Panel for Inputs ---
    sidebarPanel(
      # ... [Inputs remain the same] ... 
      sliderInput("theta_val", "\\(\\theta\\) (True parameter value)", min = -3, max = 3, value = 0, step = 0.1),
      sliderInput("x_val", "\\(x\\) (Observed data)", min = -3, max = 3, value = 1.5, step = 0.1), # 
      numericInput("prior_mean", "Prior mean (\\(\\mu_0\\))", value = 0, step = 0.1),
      numericInput("prior_sd", "Prior SD (\\(\\sigma_0\\))", value = 1, min = 0.1, step = 0.1),
      numericInput("like_sd", "Likelihood SD (\\(\\sigma\\))", value = 1, min = 0.1, step = 0.1),
      numericInput("theta_threshold", "Region Threshold \\(\\theta_0\\) (for P(\\(\\theta > \\theta_0\\)))", value = 0, step = 0.1)
    ),

    # --- Main Panel for Outputs (with Tabs) ---
    mainPanel(
      tabsetPanel(
        # --- Tab 1: Visualizations and Results ---
        tabPanel("視覺化與結果", # 
                 fluidRow(
                   column(6,
                          h5("機率分佈 (Fixed \\(\\theta\\))"), # 
                          # Explanation for Probability Plot
                          p("此圖顯示：給定一個'真實'的參數值 \\(\\theta\\)（由滑桿設定），觀察到不同數據值 \\(x\\) 的機率密度 \\( P(x \\mid \\theta) \\)。紅線標示出實際觀察到的數據點 \\(x\\)（由滑桿設定）。"), # 
                          plotOutput("prob_plot") # 
                   ),
                   column(6,
                          h5("似然函數 (Fixed \\(x\\))"), # 
                          # Explanation for Likelihood Plot
                          p("此圖顯示：給定我們觀察到的數據 \\(x\\)（由滑桿設定），不同參數值 \\(\\theta\\) 的「似然性」 \\( L(\\theta \\mid x) \\)。它衡量了不同的 \\(\\theta\\) 值對觀察數據的解釋程度。藍線標示出假設的'真實' \\(\\theta\\) 值（由滑桿設定）。"), # 
                          plotOutput("like_plot") # 
                   )
                 ),
                 fluidRow(
                   column(12, # Adjusted column width to full 
                          h5("後驗分佈"), # 
                          # Explanation for Posterior Plot
                          p("此圖整合了先驗信念（灰色虛線）和數據似然性（橘色點線，已縮放以便比較）來產生後驗分布（深綠色實線），代表了觀察數據後我們對 \\(\\theta\\) 的更新信念。圖中也顯示了："), # 
                          tags$ul(
                            tags$li("先驗均值（灰色垂直虛線）與後驗均值（深綠色垂直實線）。"), # 
                            tags$li("95% 置信區間或最高後驗密度區間（HPDI，綠色陰影區域）：根據後驗分布，我們有95%的信心認為真實 \\(\\theta\\) 在此區間內。"), # 
                            tags$li("傳統的95%信賴區間（紅色垂直虛線）：基於數據 \\(x\\) 和概似標準差 \\(\\sigma\\) 計算得出。"), # 
                            tags$li("後驗機率 \\( P(\\theta > \\theta_0) \\) 的區域（紅色陰影區域），其中 \\(\\theta_0\\) 是由使用者設定的閾值（紫色垂直點線）。圖標題顯示了計算出的機率值。") # 
                          ),
                          plotOutput("post_plot"), # 
                          # Explanation for Summary Output
                          p("以下為後驗分析的具體數值結果："),
                          uiOutput("post_summary") # 
                   )
                 )
        ),
        # --- Tab 2: Conceptual Explanation ---
        tabPanel("概念解釋", # 
                 h4("貝氏推論：知識的更新"), # 
                 p("這個應用程式展示了貝氏統計和頻率論統計在推論上的異同。貝氏統計的核心思想是利用觀察到的數據來更新我們對未知參數 \\(\\theta\\) 的信念或知識。"), # 

                 h5("核心概念"), # 
                 p(strong("機率"), "：在給定一個", strong("假設"), "的參數值 \\(\\theta\\) 的情況下，觀察到特定數據 \\(x\\) 的可能性有多大？這寫作 \\( P(x \\mid \\theta) \\)。它關注的是從一個", strong("假定的"), "模型 (\\(\\theta\\)) 來預測數據 (\\(x\\))。"), # 
                 p(strong("似然"), "：當我們已經觀察到數據 \\(x\\) 時，哪個（或哪些）參數值 \\(\\theta\\) 最能", strong("解釋"), "或", strong("符合"), "這個數據？這寫作 \\( L(\\theta \\mid x) \\)，它與 \\( P(x \\mid \\theta) \\) 在數學形式上相關，但觀點不同：數據是固定的，我們評估不同的 \\(\\theta\\)。它衡量的是不同 \\(\\theta\\) 值對觀察數據的", strong("合理性"), "。"), # 
                 p(strong("貝氏視角"), "：在貝氏框架下，\\(\\theta\\) 被視為一個我們不確定的未知量。我們對 \\(\\theta\\) 的", strong("信念或知識程度"), "可以用一個機率分布來表示（這就是", strong("先驗分布"), " \\( P(\\theta) \\) ）。當我們獲得數據 \\(x\\) 後，我們利用貝氏定理來", strong("更新"), "我們對 \\(\\theta\\) 的信念，得到", strong("後驗分布"), "（\\( P(\\theta \\mid x) \\)）："), # 
                 p(style = "margin-left:20px;", "\\( P(\\theta \\mid x) \\propto L(\\theta \\mid x) \\cdot P(\\theta) \\)"), # 
                 p(style = "margin-left:20px;", "(後驗信念 正比於 數據似然性 乘以 先驗信念)"), # 

                 h5("概念總結"), # 
                 tags$ul( # 
                   tags$li(strong("機率"), ": 從一個", strong("假設"), "的 \\(\\theta\\) 出發，預測數據 \\(x\\) 的", strong("可能性"), "。"), # 
                   tags$li(strong("似然"), ": 從", strong("觀察"), "到的數據 \\(x\\) 出發，評估哪個 \\(\\theta\\) 更能", strong("解釋"), "這個數據。"), # 
                   tags$li(strong("後驗"), ": 結合我們", strong("原本的信念"), "和", strong("數據提供的似然率"), "，得到", strong("更新後的信念"), "。") # 
                 ),

                 h5("類比：醫生診斷"), # 
                 tags$ul( # 
                    tags$li(strong("先驗"), ": 根據病人的基本情況和普遍醫學知識，醫生對病人可能患有哪些疾病有一個初步的判斷（例如，感冒的可能性 60%，流感 30%，其他 10%）。"), # 
                    tags$li(strong("似然"), ": 病人做了檢查（數據 \\(x\\)，例如發燒 39 度）。醫生會評估：如果病人得的是感冒，出現這種高燒的可能性有多大？如果得的是流感，可能性又有多大？"), # 
                    tags$li(strong("後驗"), ": 結合初步判斷（先驗）和檢查結果（似然），醫生更新診斷，可能認為流感的可能性現在更高了（例如，感冒 20%，流感 70%，其他 10%）。") # 
                 ),

                 h5("頻率派 vs. 貝氏統計：本質差異"), # 
                 tags$ul( # 
                   tags$li(strong("頻率派"), ": 參數 \\(\\theta\\) 是固定未知的常數，重點是使用重複抽樣的方法評估推論的可靠性，例如信賴區間與假設檢定。"), # 
                   tags$li(strong("貝氏統計"), ": 參數 \\(\\theta\\) 是不確定的，我們以機率分布表示我們對 \\(\\theta\\) 的信念，並透過貝氏定理更新這種信念。") # 
                 ),

                 h5("貝氏定理"), # 
                 p("貝氏定理提供了一種結合我們原有知識與新數據的方式："), # 
                 p("$$ P(\\theta \\mid x) = \\frac{P(x \\mid \\theta) P(\\theta)}{P(x)} $$"), # 
                 tags$ul( # 
                   tags$li("\\( P(\\theta) \\)：先驗分布，表示在觀察數據前我們對 \\(\\theta\\) 的信念"), # 
                   tags$li("\\( P(x \\mid \\theta) \\)：似然函數，表示在不同 \\(\\theta\\) 下觀察到數據的機率"), # 
                   tags$li("\\( P(x) \\)：邊際似然，對所有可能 \\(\\theta\\) 加總後觀察到 \\(x\\) 的總機率，確保後驗為有效機率分布"), # 
                   tags$li("\\( P(\\theta \\mid x) \\)：後驗分布，表示觀察數據後，我們對 \\(\\theta\\) 的更新信念") # 
                 ),

                 h5("貝氏因子"), # 
                 p("貝氏因子用來比較兩個模型（或假設）對觀察到的資料的支持程度，定義為："), # 
                 p("$$ BF_{10} = \\frac{P(x \\mid H_1)}{P(x \\mid H_0)} $$"), # 
                 tags$ul( # 
                   tags$li("\\( H_1 \\)：替代假說，例如 \\(\\theta > 0\\)"), # 
                   tags$li("\\( H_0 \\)：零假說，例如 \\(\\theta = 0\\)"), # 
                   tags$li("若 \\( BF_{10} > 1 \\)，表示數據對 \\( H_1 \\) 的支持比 \\( H_0 \\) 強 (BF10 =3 弱, BF10 =10 強), 越大代表越強支持；若 \\( BF_{10} < 1 (BF10 =1/3 弱, BF10 =1/10 強) \\)，則支持 \\( H_0 \\)。") # 
                 ),
                 p("這個指標是貝氏推論用來量化證據強度的重要工具，和頻率派的 p 值在理念上不同。"), # 

                 h5("表格比較"), # 
                 tags$table(border = '1', style = 'border-collapse: collapse; width: 100%; text-align: center;', # 
                    tags$thead( # 
                      tags$tr( # 
                        tags$th("概念"), tags$th("固定什麼？"), tags$th("變化/評估什麼？"), tags$th("代表什麼？") # 
                      )
                    ),
                    tags$tbody( # 
                      tags$tr( # 
                        tags$td("機率"), tags$td("\\(\\theta\\) (假設的模型/參數)"), tags$td("\\(x\\) (潛在的數據)"), tags$td("在特定假設下，預測數據出現的可能性") # 
                      ),
                      tags$tr( # 
                        tags$td("似然"), tags$td("\\(x\\) (觀察到的數據)"), tags$td("\\(\\theta\\) (不同的參數值)"), tags$td("評估不同參數值對觀察數據的解釋力/合理性") # 
                      ),
                      tags$tr( # 
                        tags$td("後驗"), tags$td("\\(x\\) (觀察到的數據)"), tags$td("\\(\\theta\\) (以分布形式表示的信念)"), tags$td("經過數據更新後，我們對 \\(\\theta\\) 的最終信念狀態") # 
                      )
                    )
                 )
        ) # End tabPanel Conceptual Explanation 
      ) # End tabsetPanel 
    ) # End mainPanel 
  ) # End sidebarLayout 
) # End fluidPage 

# --- Server Logic ---
server <- function(input, output, session) { # 

  # --- Probability Plot ---
  output$prob_plot <- renderPlot({ # 
    # ... [Probability Plot code remains the same] ... 
    theta <- input$theta_val
    like_sd <- input$like_sd
    x_vals <- seq(theta - 4*like_sd, theta + 4*like_sd, length.out = 400)
    y_vals <- dnorm(x_vals, mean = theta, sd = like_sd)

    plot(x_vals, y_vals, type = "l", lwd = 2, col = "steelblue", # 
         main = paste0("P(x | θ = ", theta, ", σ = ", like_sd, ")"), # 
         ylab = "Density", xlab = "x", ylim = c(0, max(y_vals) * 1.1))
    abline(v = input$x_val, col = "red", lty = 2) # 
    text(input$x_val, max(y_vals) * 0.1, paste0("Observed x = ", input$x_val), pos = ifelse(input$x_val > theta, 2, 4), col = "red") # 
  })

  # --- Likelihood Plot ---
  output$like_plot <- renderPlot({ # 
    # ... [Likelihood Plot code remains the same] ... 
    x <- input$x_val
    like_sd <- input$like_sd
    # Adjust range based on prior and likelihood SDs for better visualization 
    prior_mean <- input$prior_mean # 
    prior_sd <- input$prior_sd # 
    plot_range_sd <- max(prior_sd, like_sd) * 5 # Heuristic range 
    theta_vals <- seq(min(x - plot_range_sd, prior_mean - plot_range_sd), # 
                      max(x + plot_range_sd, prior_mean + plot_range_sd), length.out = 400)
    like_vals <- dnorm(x, mean = theta_vals, sd = like_sd) # 

    plot(theta_vals, like_vals, type = "l", lwd = 2, col = "darkorange", # 
         main = paste0("L(θ | x = ", x, ", σ = ", like_sd, ")"), # 
         ylab = "Likelihood (Unnormalized)", xlab = "θ", ylim = c(0, max(like_vals) * 1.1))
    abline(v = input$theta_val, col = "blue", lty = 2) # 
    text(input$theta_val, max(like_vals) * 0.9, paste0("True θ = ", input$theta_val), pos = 4, col = "blue") # 
  })

  # Reactive expression for posterior parameters
  posterior_params <- reactive({ # 
    # ... [Posterior Parameters calculation remains the same] ... 
    x <- input$x_val
    sigma <- input$like_sd
    mu0 <- input$prior_mean
    sigma0 <- input$prior_sd # 

    prec_prior <- 1 / sigma0^2 # 
    prec_like <- 1 / sigma^2 # 
    prec_post <- prec_prior + prec_like # 
    mu_post <- (mu0 * prec_prior + x * prec_like) / prec_post # 
    sigma_post <- sqrt(1 / prec_post) # 

    list(mu = mu_post, sd = sigma_post, prec = prec_post) # 
  })

  # --- Posterior Plot ---
  output$post_plot <- renderPlot({ # 
    # ... [Posterior Plot code remains the same] ... 
    params <- posterior_params()
    mu_post <- params$mu
    sigma_post <- params$sd
    theta0 <- input$theta_threshold

    x <- input$x_val # 
    sigma <- input$like_sd # 
    mu0 <- input$prior_mean # 
    sigma0 <- input$prior_sd # 

    # Determine plot range dynamically 
    plot_min <- min(mu_post - 4*sigma_post, mu0 - 4*sigma0, x - 4*sigma) # 
    plot_max <- max(mu_post + 4*sigma_post, mu0 + 4*sigma0, x + 4*sigma) # 
    theta_vals <- seq(plot_min, plot_max, length.out = 400) # 

    # Calculate densities 
    prior_vals <- dnorm(theta_vals, mu0, sigma0) # 
    # Scale likelihood to be visible with prior/posterior 
    like_vals_raw <- dnorm(x, mean = theta_vals, sd = sigma) # 
    # Simple scaling: scale max likelihood to match max prior or posterior density 
    max_dens <- max(c(prior_vals, dnorm(theta_vals, mu_post, sigma_post))) # 
    like_vals_scaled <- if (max(like_vals_raw) > 1e-6) { # Avoid division by zero or near-zero 
                       like_vals_raw / max(like_vals_raw) * max_dens # 
                     } else {
                       rep(0, length(theta_vals)) # Handle case where likelihood is flat zero 
                     }
    post_vals <- dnorm(theta_vals, mu_post, sigma_post) # 

    # Calculate HPD interval and Frequentist CI 
    post_sample <- rnorm(10000, mu_post, sigma_post) # Sample for HPD 
    hpd_ci <- hdi(post_sample, credMass = 0.95) # 
    freq_lower <- x - 1.96 * sigma # 
    freq_upper <- x + 1.96 * sigma # 

    # Calculate posterior probability P(theta > theta0) 
    prob_gt_theta0 <- 1 - pnorm(theta0, mean = mu_post, sd = sigma_post) # 

    # Create the plot 
    plot(theta_vals, post_vals, type = "l", lwd = 2, col = "darkgreen", # 
         ylab = "Density / Scaled Likelihood", xlab = "θ",
         main = paste0("Posterior Distribution: P(θ > ", theta0, ") ≈ ", round(prob_gt_theta0, 3)), # 
         ylim = c(0, max(c(prior_vals, like_vals_scaled, post_vals)) * 1.1)) # Adjust y-axis limit 

    # Add shaded regions and lines 
    # Shaded region for P(theta > theta0) 
    region_gt_indices <- which(theta_vals >= theta0) # 
    if (length(region_gt_indices) > 0) { # 
      polygon(c(theta_vals[region_gt_indices], rev(theta_vals[region_gt_indices])), # 
              c(post_vals[region_gt_indices], rep(0, length(region_gt_indices))), # 
              col = rgb(1, 0, 0, 0.2), border = NA) # 
    }
    # Shaded region for 95% HPD 
    region_hpd_indices <- which(theta_vals >= hpd_ci[1] & theta_vals <= hpd_ci[2]) # 
     if (length(region_hpd_indices) > 0) { # 
       polygon(c(theta_vals[region_hpd_indices], rev(theta_vals[region_hpd_indices])), # 
               c(post_vals[region_hpd_indices], rep(0, length(region_hpd_indices))), # 
               col = rgb(0, 1, 0, 0.2), border = NA) # 
     }

    # Add lines for distributions and CIs 
    lines(theta_vals, prior_vals, col = "gray", lwd = 2, lty = 2) # 
    lines(theta_vals, like_vals_scaled, col = "darkorange", lty = 3, lwd = 2) # 
    lines(theta_vals, post_vals, col = "darkgreen", lwd = 2) # Redraw posterior line on top 

    # Add vertical lines for means, CIs, threshold 
    abline(v = mu_post, col = "darkgreen", lty = 1, lwd=1) # Posterior mean 
    abline(v = mu0, col = "gray", lty = 2)                 # Prior mean 
    abline(v = freq_lower, col = "red", lty = 2)           # Freq CI lower 
    abline(v = freq_upper, col = "red", lty = 2)           # Freq CI upper 
    abline(v = theta0, col = "purple", lty = 4)            # Threshold 

    # Add legend 
    legend("topright", # 
           legend = c("Posterior", "Prior", "Likelihood (scaled)", "95% HPD Region", paste0("Region θ > ", theta0), "Posterior Mean", "Prior Mean", "Freq. 95% CI", "Threshold θ₀"), # 
           col = c("darkgreen", "gray", "darkorange", rgb(0,1,0,0.4), rgb(1,0,0,0.4), "darkgreen", "gray", "red", "purple"), # 
           lty = c(1, 2, 3, NA, NA, 1, 2, 2, 4), # 
           lwd = c(2, 2, 2, NA, NA, 1, 1, 1, 1), # 
           fill = c(NA, NA, NA, rgb(0,1,0,0.2), rgb(1,0,0,0.2), NA, NA, NA, NA), # 
           border = NA, bty = "n", cex=0.9) # Smaller font size for legend 
  })

  # --- Posterior Summary Text ---
  output$post_summary <- renderUI({ # 
      # ... [Posterior Summary calculation remains the same] ... 
      params <- posterior_params()
      mu_post <- params$mu
      sigma_post <- params$sd
      theta0 <- input$theta_threshold
      x <- input$x_val
      sigma <- input$like_sd

      # Calculate intervals and probability again (or retrieve from reactive) 
      post_sample <- rnorm(10000, mu_post, sigma_post) # 
      hpd_ci <- hdi(post_sample, credMass = 0.95) # 
      freq_ci <- c(x - 1.96 * sigma, x + 1.96 * sigma) # 
      prob_gt_theta0 <- 1 - pnorm(theta0, mu_post, sigma_post) # 

      # Format output as HTML 
      HTML(sprintf( # 
         "<b>Posterior Mean:</b> %.3f<br>
          <b>Posterior SD:</b> %.3f<br>
          <b>95%% Highest Posterior Density Interval (HPDI):</b> [%.3f, %.3f]<br>
          <b>Frequentist 95%% Confidence Interval (CI):</b> [%.3f, %.3f]<br>
          <b>Posterior Probability P(θ > %.2f):</b> %.3f", # 
         mu_post, sigma_post, hpd_ci[1], hpd_ci[2], # 
         freq_ci[1], freq_ci[2], theta0, prob_gt_theta0 # 
      ))
  })

} # End server 

# --- Run the Application ---
shinyApp(ui = ui, server = server) # 
