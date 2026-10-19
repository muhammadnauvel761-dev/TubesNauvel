# ==============================================================================
# 0. OTOMATIS INSTAL & PANGGIL PACKAGE
# ==============================================================================
paket_dibutuhkan <- c(
  "sf", "terra", "raster", "dplyr", 
  "shiny", "shinydashboard", "leaflet", "leaflet.extras", "DT", "plotly"
)

paket_baru <- paket_dibutuhkan[!(paket_dibutuhkan %in% installed.packages()[,"Package"])]
if(length(paket_baru)) install.packages(paket_baru)

library(sf)
library(terra)
library(raster)
library(dplyr)
library(shiny)
library(shinydashboard)
library(leaflet)
library(leaflet.extras)
library(DT)
library(plotly)

# ==============================================================================
# 1. PENGATURAN PROJ & LINGKUNGAN KERJA
# ==============================================================================
Sys.unsetenv("PROJ_LIB")
Sys.unsetenv("PROJ_DATA")
Sys.unsetenv("GDAL_DATA")

# ==============================================================================
# 2. MENENTUKAN FOLDER & JALUR DATA
# ==============================================================================
folder_tubes <- "D:/Kuliah SMT5/Pemspa P Neo/TUBES"
setwd(folder_tubes)

folder_citra <- file.path(
  folder_tubes,
  "S2C_MSIL2A_20260529T022531_N0512_R046_T49MFM_20260529T061812.SAFE/S2C_MSIL2A_20260529T022531_N0512_R046_T49MFM_20260529T061812.SAFE/GRANULE/L2A_T49MFM_A009025_20260529T024756/IMG_DATA"
)

folder_batas <- file.path(folder_tubes, "batas_malang_fiks")
folder_kecamatan <- file.path(folder_tubes, "batas_kecamatan_malang")

# Memeriksa folder
stopifnot(dir.exists(folder_citra))
stopifnot(dir.exists(folder_batas))

# ==============================================================================
# 3. PEMROSESAN CITRA SATELIT & CROPPING
# ==============================================================================
files <- list.files(
  folder_citra,
  pattern = "\\.jp2$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

ambil_band <- function(kode) {
  hasil <- grep(kode, files, value = TRUE)
  if (length(hasil) != 1) {
    stop(paste("Jumlah file untuk band", kode, "=", length(hasil)))
  }
  rast(hasil)
}

B03 <- ambil_band("R10m/T49MFM_20260529T022531_B03_10m")
B04 <- ambil_band("R10m/T49MFM_20260529T022531_B04_10m")
B08 <- ambil_band("R10m/T49MFM_20260529T022531_B08_10m")
B11 <- ambil_band("R20m/T49MFM_20260529T022531_B11_20m")

# Batas Malang
shp_files <- list.files(folder_batas, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
malang <- st_read(shp_files[1]) %>% st_make_valid()
batas <- terra::vect(malang)

if (!same.crs(batas, B11)) {
  batas <- terra::project(batas, terra::crs(B11))
}

# Samakan resolusi band & penggabungan
B03_new <- resample(B03, B11, method = "bilinear")
B04_new <- resample(B04, B11, method = "bilinear")
B08_new <- resample(B08, B11, method = "bilinear")

citra <- c(B03_new, B04_new, B08_new, B11)
names(citra) <- c("B03", "B04", "B08", "B11")

# Cropping & Masking
citra_crop <- crop(citra, batas)
citra_mask <- mask(citra_crop, batas)

# Buat folder output
folder_output <- file.path(folder_tubes, "output")
dir.create(folder_output, showWarnings = FALSE, recursive = TRUE)

writeRaster(citra_mask, file.path(folder_output, "citra_malang.tif"), overwrite = TRUE)

# ==============================================================================
# 4. PERHITUNGAN INDEKS SPEKTRAL (NDVI, NDBI, NDWI, NDCI)
# ==============================================================================
NDVI <- terra::ifel((citra_mask$B08 + citra_mask$B04) == 0, NA, (citra_mask$B08 - citra_mask$B04) / (citra_mask$B08 + citra_mask$B04))
names(NDVI) <- "NDVI"
writeRaster(NDVI, file.path(folder_output, "NDVI_malang.tif"), overwrite = TRUE)

NDBI <- terra::ifel((citra_mask$B11 + citra_mask$B08) == 0, NA, (citra_mask$B11 - citra_mask$B08) / (citra_mask$B11 + citra_mask$B08))
names(NDBI) <- "NDBI"
writeRaster(NDBI, file.path(folder_output, "NDBI_malang.tif"), overwrite = TRUE)

NDWI <- terra::ifel((citra_mask$B03 + citra_mask$B08) == 0, NA, (citra_mask$B03 - citra_mask$B08) / (citra_mask$B03 + citra_mask$B08))
names(NDWI) <- "NDWI"
writeRaster(NDWI, file.path(folder_output, "NDWI_malang.tif"), overwrite = TRUE)

Red <- citra_mask$B04
NIR <- citra_mask$B08
NDCI <- terra::ifel((Red + NIR) == 0, NA, (Red - NIR) / (Red + NIR))
names(NDCI) <- "NDCI"
writeRaster(NDCI, file.path(folder_output, "NDCI_Cloud_Malang.tif"), overwrite = TRUE)

# ==============================================================================
# 5. KLASIFIKASI TUTUPAN LAHAN 4 KELAS
# ==============================================================================
klasifikasi <- terra::ifel(NDCI > 0.05, 1,
                           terra::ifel(NDWI > 0, 2,
                                       terra::ifel(NDVI > 0.3, 3,
                                                   terra::ifel(NDBI > 0, 4, 5)
                                       )
                           )
)

levels(klasifikasi) <- data.frame(
  ID = c(1, 2, 3, 4, 5),
  Kategori = c("Awan", "Perairan", "Vegetasi", "Bangunan", "Lainnya")
)

writeRaster(klasifikasi, file.path(folder_output, "Klasifikasi_Tutupan_Lahan.tif"), overwrite = TRUE)

# ==============================================================================
# 6. ZONAL STATISTICS PER KECAMATAN
# ==============================================================================
shp_kecamatan <- list.files(folder_kecamatan, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
kecamatan <- st_read(shp_kecamatan[1]) %>% st_make_valid()
kecamatan <- kecamatan[!st_is_empty(kecamatan), ]
kecamatan <- st_transform(kecamatan, crs = terra::crs(NDVI))
kec <- terra::vect(kecamatan)

samakan_grid <- function(r, acuan) {
  if (!terra::same.crs(r, acuan)) {
    r <- terra::project(r, acuan, method = "bilinear")
  } else if (!terra::compareGeom(r, acuan, stopOnError = FALSE)) {
    r <- terra::resample(r, acuan, method = "bilinear")
  }
  r
}

NDBI_fix <- samakan_grid(NDBI, NDVI)
NDWI_fix <- samakan_grid(NDWI, NDVI)
NDCI_fix <- samakan_grid(NDCI, NDVI)

indeks_zonal <- c(NDVI, NDBI_fix, NDWI_fix, NDCI_fix)
names(indeks_zonal) <- c("NDVI", "NDBI", "NDWI", "NDCI")

kolom_kec <- if("WADMKC" %in% names(kecamatan)) "WADMKC" else names(kecamatan)[1]

zonal <- data.frame(Kecamatan = kecamatan[[kolom_kec]])
statistik <- list(Mean = mean, Minimum = min, Maximum = max, Median = median, SD = sd)

for (nama in names(indeks_zonal)) {
  for (stat in names(statistik)) {
    hasil <- terra::extract(indeks_zonal[[nama]], kec, fun = statistik[[stat]], na.rm = TRUE)
    zonal[[paste(nama, stat, sep = "_")]] <- hasil[[2]]
  }
}

write.csv(zonal, file.path(folder_output, "Zonal_Statistics_Malang.csv"), row.names = FALSE)

kecamatan_merged <- kecamatan %>%
  left_join(zonal, by = setNames("Kecamatan", kolom_kec))

# ==============================================================================
# 7. WEB GIS DASHBOARD INTERAKTIF (R SHINY & LEAFLET)
# ==============================================================================
kec_leaflet <- st_transform(kecamatan_merged, crs = 4326)

ndvi_leaf <- raster(NDVI)
ndbi_leaf <- raster(NDBI_fix)
ndwi_leaf <- raster(NDWI_fix)
klas_leaf <- raster(klasifikasi)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "WebGIS Indeks Spektral Malang"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Peta Interaktif", tabName = "m_peta", icon = icon("map")),
      menuItem("Statistik Zonal", tabName = "m_statistik", icon = icon("chart-bar")),
      menuItem("Tabel Zonal Statistics", tabName = "m_tabel", icon = icon("table"))
    )
  ),
  
  dashboardBody(
    tags$style(type = 'text/css', '#map_indeks {height: calc(100vh - 80px) !important;}'),
    
    tabItems(
      # TAB 1: PETA INTERAKTIF
      tabItem(
        tabName = "m_peta",
        leafletOutput('map_indeks', height = 700)
      ),
      
      # TAB 2: STATISTIK ZONAL PER KECAMATAN (DESAIN SESUAI RUJUKAN)
      tabItem(
        tabName = "m_statistik",
        fluidRow(
          column(width = 12, h3("Statistik Zonal per Kecamatan", style = "margin-top:0px; font-weight:bold;")),
          column(width = 12, p("Pilih indeks, lihat peringkatnya, lalu telusuri angka rincinya.", style = "color: #666;"))
        ),
        
        # 4 KARTU METRIK (VALUE BOXES)
        fluidRow(
          valueBoxOutput("box_rata", width = 3),
          valueBoxOutput("box_max", width = 3),
          valueBoxOutput("box_min", width = 3),
          valueBoxOutput("box_total", width = 3)
        ),
        
        br(),
        
        # PANEL KONTROL KIRI & DIAGRAM KANAN
        fluidRow(
          box(
            title = "Pilihan Tampilan", status = "primary", solidHeader = TRUE, width = 4,
            selectInput("pilih_indeks", "Indeks Spektral:",
                        choices = c("NDVI - Vegetasi" = "NDVI_Mean",
                                    "NDBI - Area Terbangun" = "NDBI_Mean",
                                    "NDWI - Perairan" = "NDWI_Mean",
                                    "NDCI - Awan" = "NDCI_Mean"),
                        selected = "NDVI_Mean"),
            uiOutput("slider_filter_ui")
          ),
          
          box(
            title = "Peringkat Rata-rata tiap Kecamatan", status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("plot_peringkat", height = "500px")
          )
        )
      ),
      
      # TAB 3: TABEL DATA
      tabItem(
        tabName = "m_tabel",
        fluidRow(
          box(
            title = "Tabel Zonal Statistics Kabupaten Malang", 
            status = "primary", solidHeader = TRUE, width = 12,
            DTOutput('tabel_zonal')
          )
        )
      )
    )
  )
)

server <- function(input, output, session){
  
  # --- UI DYNAMIC SLIDER BASED ON INDEX VALUE ---
  output$slider_filter_ui <- renderUI({
    req(input$pilih_indeks)
    val_min <- floor(min(zonal[[input$pilih_indeks]], na.rm = TRUE) * 100) / 100
    val_max <- ceiling(max(zonal[[input$pilih_indeks]], na.rm = TRUE) * 100) / 100
    
    sliderInput("filter_val", "Saring nilai rata-rata:",
                min = val_min, max = val_max,
                value = c(val_min, val_max), step = 0.01)
  })
  
  # --- FILTERED DATASET ---
  data_filtered <- reactive({
    req(input$pilih_indeks, input$filter_val)
    zonal %>%
      filter(.data[[input$pilih_indeks]] >= input$filter_val[1] & 
               .data[[input$pilih_indeks]] <= input$filter_val[2])
  })
  
  # --- 1. VALUE BOX: RATA-RATA KESELURUHAN ---
  output$box_rata <- renderValueBox({
    req(input$pilih_indeks)
    val <- mean(zonal[[input$pilih_indeks]], na.rm = TRUE)
    valueBox(
      value = sprintf("%.3f", val),
      subtitle = "Rata-rata Keseluruhan",
      icon = icon("chart-line"),
      color = "light-blue"
    )
  })
  
  # --- 2. VALUE BOX: TERTINGGI ---
  output$box_max <- renderValueBox({
    req(input$pilih_indeks)
    df <- zonal %>% arrange(desc(.data[[input$pilih_indeks]])) %>% slice(1)
    valueBox(
      value = sprintf("%.3f", df[[input$pilih_indeks]]),
      subtitle = paste("Tertinggi -", df$Kecamatan),
      icon = icon("arrow-up"),
      color = "green"
    )
  })
  
  # --- 3. VALUE BOX: TERENDAH ---
  output$box_min <- renderValueBox({
    req(input$pilih_indeks)
    df <- zonal %>% arrange(.data[[input$pilih_indeks]]) %>% slice(1)
    valueBox(
      value = sprintf("%.3f", df[[input$pilih_indeks]]),
      subtitle = paste("Terendah -", df$Kecamatan),
      icon = icon("arrow-down"),
      color = "orange"
    )
  })
  
  # --- 4. VALUE BOX: JUMLAH KECAMATAN ---
  output$box_total <- renderValueBox({
    valueBox(
      value = nrow(zonal),
      subtitle = "Kecamatan dianalisis",
      icon = icon("map-marked-alt"),
      color = "purple"
    )
  })
  
  # --- DIAGRAM BATANG PERINGKAT HORISONTAL (PLOTLY) ---
  output$plot_peringkat <- renderPlotly({
    df <- data_filtered()
    req(nrow(df) > 0)
    
    col_name <- input$pilih_indeks
    df_sorted <- df %>% arrange(.data[[col_name]]) # Terurut untuk horizontal bar
    
    # Menentukan warna batang berdasarkan jenis indeks
    bar_color <- switch(col_name,
                        "NDVI_Mean" = "#2ecc71",
                        "NDBI_Mean" = "#e67e22",
                        "NDWI_Mean" = "#3498db",
                        "NDCI_Mean" = "#95a5a6")
    
    plot_ly(
      data = df_sorted,
      x = ~get(col_name),
      y = ~reorder(Kecamatan, get(col_name)),
      type = 'bar',
      orientation = 'h',
      marker = list(color = bar_color),
      text = ~sprintf("%.3f", get(col_name)),
      textposition = 'outside'
    ) %>%
      layout(
        xaxis = list(title = "Rata-rata"),
        yaxis = list(title = ""),
        margin = list(l = 120, r = 50, t = 20, b = 40)
      )
  })
  
  # --- PETA INTERAKTIF LAYER LEAFLET ---
  output$map_indeks <- renderLeaflet({
    pal_ndvi <- colorNumeric("YlGn", values(ndvi_leaf), na.color = "transparent")
    pal_ndbi <- colorNumeric("YlOrRd", values(ndbi_leaf), na.color = "transparent")
    pal_ndwi <- colorNumeric("Blues", values(ndwi_leaf), na.color = "transparent")
    
    popup_kec <- paste0(
      "<strong>Kecamatan: </strong>", kec_leaflet[[kolom_kec]], "<br>",
      "<hr style='margin: 5px 0;'>",
      "<b>Rata-rata NDVI (Vegetasi): </b>", round(kec_leaflet$NDVI_Mean, 3), "<br>",
      "<b>Rata-rata NDBI (Bangunan): </b>", round(kec_leaflet$NDBI_Mean, 3), "<br>",
      "<b>Rata-rata NDWI (Air): </b>", round(kec_leaflet$NDWI_Mean, 3), "<br>",
      "<b>Rata-rata NDCI (Awan): </b>", round(kec_leaflet$NDCI_Mean, 3)
    )
    
    leaflet() %>%
      addTiles(group = "OpenStreetMap") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
      addRasterImage(ndvi_leaf, colors = pal_ndvi, opacity = 0.8, group = "NDVI (Vegetasi)") %>%
      addRasterImage(ndbi_leaf, colors = pal_ndbi, opacity = 0.8, group = "NDBI (Bangunan)") %>%
      addRasterImage(ndwi_leaf, colors = pal_ndwi, opacity = 0.8, group = "NDWI (Air)") %>%
      addPolygons(
        data = kec_leaflet,
        fillColor = "transparent",
        weight = 1.5,
        color = "white",
        popup = popup_kec,
        group = "Batas Kecamatan"
      ) %>%
      addLayersControl(
        baseGroups = c("OpenStreetMap", "Esri World Imagery"),
        overlayGroups = c("NDVI (Vegetasi)", "NDBI (Bangunan)", "NDWI (Air)", "Batas Kecamatan"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  output$tabel_zonal <- renderDT({
    datatable(zonal, options = list(pageLength = 10, scrollX = TRUE))
  })
}

# ==============================================================================
# 8. MENJALANKAN APLIKASI WEB GIS SHINY
# ==============================================================================
shinyApp(ui, server)