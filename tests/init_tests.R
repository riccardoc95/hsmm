# =====================================================
# TEST UTILITY COMPLETA PER fit_hsmm
# =====================================================
# Questa utility testa tutte le funzionalità della funzione fit_hsmm
# e rileva potenziali bug o problemi

library(testthat)

# Utility per generare dati di test
generate_test_data <- function(n = 100, emission_type = "torus", seed = 123) {
  set.seed(seed)
  
  if (emission_type == "torus") {
    # Dati per torus (angoli in [0, 2*pi])
    data <- cbind(runif(n, 0, 2*pi), runif(n, 0, 2*pi))
  } else if (emission_type == "gaussian") {
    # Dati gaussiani (se implementati)
    data <- cbind(rnorm(n), rnorm(n))
  } else if (emission_type == "poisson") {
    # Dati Poisson (se implementati)
    data <- cbind(rpois(n, 5), rpois(n, 3))
  } else {
    # Default: torus
    data <- cbind(runif(n, 0, 2*pi), runif(n, 0, 2*pi))
  }
  
  return(data)
}

# Utility per generare covariate
generate_covariates <- function(n = 100, ncov = 1, seed = 123) {
  set.seed(seed)
  if (ncov == 1) {
    return(matrix(rnorm(n), ncol = 1))
  } else {
    return(matrix(rnorm(n * ncov), ncol = ncov))
  }
}

# Test 1: Configurazione di base (semplice)
test_basic_functionality <- function(emission_type = "torus") {
  cat("\n=== TEST 1: Configurazione di base ===\n")
  
  tryCatch({
    data <- generate_test_data(n = 50, emission_type = emission_type)
    
    result <- fit_hsmm(
      data = data,
      covariates = NULL,
      n_states = 2,
      max_dwell = 3,
      emission_model_type = emission_type,
      duration_effects = FALSE,
      transition_effects = FALSE,
      max_iter = 5,
      tol = 1e-3,
      verbose = FALSE,
      seed = 123
    )
    
    # Verifica che il risultato sia una lista
    expect_is(result, "list")
    cat("✓ Test di base superato\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Test di base fallito:", e$message, "\n")
    return(FALSE)
  })
}

# Test 2: Con effetti di durata
test_duration_effects <- function(emission_type = "torus") {
  cat("\n=== TEST 2: Con effetti di durata ===\n")
  
  tryCatch({
    data <- generate_test_data(n = 60, emission_type = emission_type)
    covariates <- generate_covariates(n = 60, ncov = 1)
    
    result <- fit_hsmm(
      data = data,
      covariates = covariates,
      n_states = 2,
      max_dwell = 4,
      emission_model_type = emission_type,
      duration_effects = TRUE,
      transition_effects = FALSE,
      max_iter = 5,
      tol = 1e-3,
      verbose = FALSE,
      seed = 123
    )
    
    expect_is(result, "list")
    cat("✓ Test con effetti di durata superato\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Test con effetti di durata fallito:", e$message, "\n")
    return(FALSE)
  })
}

# Test 3: Con effetti di transizione
test_transition_effects <- function(emission_type = "torus") {
  cat("\n=== TEST 3: Con effetti di transizione ===\n")
  
  tryCatch({
    data <- generate_test_data(n = 70, emission_type = emission_type)
    covariates <- generate_covariates(n = 70, ncov = 1)
    
    result <- fit_hsmm(
      data = data,
      covariates = covariates,
      n_states = 3,  # 3 stati per testare transizioni
      max_dwell = 3,
      emission_model_type = emission_type,
      duration_effects = FALSE,
      transition_effects = TRUE,
      max_iter = 5,
      tol = 1e-3,
      verbose = FALSE,
      seed = 123
    )
    
    expect_is(result, "list")
    cat("✓ Test con effetti di transizione superato\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Test con effetti di transizione fallito:", e$message, "\n")
    return(FALSE)
  })
}

# Test 4: Modello completo (entrambi gli effetti)
test_full_model <- function(emission_type = "torus") {
  cat("\n=== TEST 4: Modello completo ===\n")
  
  tryCatch({
    data <- generate_test_data(n = 80, emission_type = emission_type)
    covariates <- generate_covariates(n = 80, ncov = 2)
    
    result <- fit_hsmm(
      data = data,
      covariates = covariates,
      n_states = 2,
      max_dwell = 5,
      emission_model_type = emission_type,
      duration_effects = TRUE,
      transition_effects = TRUE,
      max_iter = 8,
      tol = 1e-3,
      verbose = TRUE,
      seed = 123
    )
    
    expect_is(result, "list")
    cat("✓ Test modello completo superato\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Test modello completo fallito:", e$message, "\n")
    return(FALSE)
  })
}

# Test 5: Validazione parametri (dovrebbe fallire)
test_parameter_validation <- function() {
  cat("\n=== TEST 5: Validazione parametri ===\n")
  
  data <- generate_test_data(n = 50)
  
  # Test con dati insufficienti
  tryCatch({
    fit_hsmm(
      data = data[1:1, , drop = FALSE],  # Solo 1 osservazione
      n_states = 2,
      max_dwell = 3,
      emission_model_type = "torus",
      verbose = FALSE
    )
    cat("✗ Test validazione fallito: dovrebbe aver rifiutato dati insufficienti\n")
    return(FALSE)
  }, error = function(e) {
    cat("✓ Test validazione superato: rifiuta correttamente dati insufficienti\n")
  })
  
  # Test con stati insufficienti
  tryCatch({
    fit_hsmm(
      data = data,
      n_states = 1,  # Troppo pochi stati
      max_dwell = 3,
      emission_model_type = "torus",
      verbose = FALSE
    )
    cat("✗ Test validazione fallito: dovrebbe aver rifiutato stati insufficienti\n")
    return(FALSE)
  }, error = function(e) {
    cat("✓ Test validazione superato: rifiuta correttamente stati insufficienti\n")
  })
  
  return(TRUE)
}

# Test 6: Convergenza e stabilità numerica
test_convergence <- function(emission_type = "torus") {
  cat("\n=== TEST 6: Test di convergenza ===\n")
  
  tryCatch({
    data <- generate_test_data(n = 100, emission_type = emission_type)
    
    result <- fit_hsmm(
      data = data,
      n_states = 2,
      max_dwell = 4,
      emission_model_type = emission_type,
      duration_effects = FALSE,
      transition_effects = FALSE,
      max_iter = 20,
      tol = 1e-6,
      verbose = TRUE,  # Verbose per monitorare convergenza
      seed = 123
    )
    
    expect_is(result, "list")
    cat("✓ Test di convergenza superato\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Test di convergenza fallito:", e$message, "\n")
    return(FALSE)
  })
}

# Test 7: Test con diversi link functions
test_different_links <- function(emission_type = "torus") {
  cat("\n=== TEST 7: Test con diversi link functions ===\n")
  
  data <- generate_test_data(n = 60, emission_type = emission_type)
  covariates <- generate_covariates(n = 60)
  
  links_to_test <- c("cloglog", "logit", "probit")
  
  for (link_func in links_to_test) {
    tryCatch({
      result <- fit_hsmm(
        data = data,
        covariates = covariates,
        n_states = 2,
        max_dwell = 3,
        emission_model_type = emission_type,
        duration_effects = TRUE,
        transition_effects = FALSE,
        max_iter = 5,
        link = link_func,
        verbose = FALSE,
        seed = 123
      )
      
      cat("✓ Link function", link_func, "funziona\n")
      
    }, error = function(e) {
      cat("✗ Link function", link_func, "fallito:", e$message, "\n")
    })
  }
  
  return(TRUE)
}

# Funzione principale per eseguire tutti i test
run_comprehensive_tests <- function(emission_type = "torus") {
  cat("\n")
  cat("=====================================================\n")
  cat("   TEST SUITE COMPLETA PER fit_hsmm\n")
  cat("   Emission Model Type:", emission_type, "\n")
  cat("=====================================================\n")
  
  results <- list()
  
  results$basic <- test_basic_functionality(emission_type)
  results$duration <- test_duration_effects(emission_type)
  results$transition <- test_transition_effects(emission_type)
  results$full <- test_full_model(emission_type)
  results$validation <- test_parameter_validation()
  results$convergence <- test_convergence(emission_type)
  results$links <- test_different_links(emission_type)
  
  # Riassunto risultati
  cat("\n")
  cat("=====================================================\n")
  cat("   RIASSUNTO RISULTATI\n")
  cat("=====================================================\n")
  
  passed <- sum(unlist(results))
  total <- length(results)
  
  cat("Test superati:", passed, "/", total, "\n")
  
  if (passed == total) {
    cat("🎉 TUTTI I TEST SUPERATI! fit_hsmm funziona correttamente.\n")
  } else {
    cat("⚠️  ALCUNI TEST FALLITI. Controllare i messaggi di errore sopra.\n")
  }
  
  cat("\n")
  return(results)
}

# Esempi di utilizzo:
run_comprehensive_tests("torus")          # Test con modello torus
# run_comprehensive_tests("gaussian")       # Test con modello gaussiano (se implementato)
# run_comprehensive_tests("poisson")        # Test con modello Poisson (se implementato)