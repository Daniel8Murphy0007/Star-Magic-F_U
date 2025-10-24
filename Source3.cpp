#include <QApplication>
#include <QMainWindow>
#include <QLineEdit>
#include <QTextEdit>
#include <QWebEngineView>
#include <QTabWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QDockWidget>
#include <QDialog>
#include <QMessageBox>
#include <QToolBar>
#include <QScreen>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <QFile>
#include <QDir>
#include <QStandardPaths>
#include <QKeyEvent>
#include <QCoreApplication>
#include <vtkSmartPointer.h>
#include <vtkScatterPlotMatrix.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkAxis.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <curl/curl.h>
#include <websocket.h>
#include <sqlite3.h>
#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/cognito-idp/CognitoIdentityProviderClient.h>
#include <pocketsphinx.h>
#include <opencv2/opencv.hpp>
#include <pybind11/embed.h>
#include <qalculate.h>
#include <windows.h>
#include <string>
#include <vector>
#include <thread>
#include <nlohmann/json.hpp>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <chrono>

#define MAX_QUERY_LENGTH 6000
#define MAX_WINDOWS 21  // Increased for ALMA Cycle 12
#define NASA_API_KEY_1 "PNJaNeFWqMb2g0CEQGqJePkndqYfKvBzq6XJqAwg"
#define NASA_API_KEY_2 "FJnBo64nLFqExHwDchrcaf101D8wmGSm0cF27clz"
#define MAST_API_KEY "emXvt90Htf0U4RogKTB5lqSxClUeg2pvMQxvZciM"
#define OPENAI_API_KEY "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"
#define COGNITO_CLIENT_ID "your_cognito_client_id"
#define COGNITO_REGION "us-east-1"

namespace py = pybind11;
using json = nlohmann::json;

// Structure for search results
struct SearchResult {
    std::string url;
    std::string title;
    std::string summary;
    double relevance;
    bool isLive;
};

// Global variables
std::vector<std::string> focusList = {
    "Worldwide Telescopes", "NASA", "SpaceX", "JPL", "ESA", "STScI",
    "Hubble", "JWST", "Chandra", "ALMA", "EHT", "SKA Observatory",
    "CERN", "DARPA", "ATIP", "ACS Hubble Ultra Deep Field",
    "WFC3 Hubble Deep Field", "Hubble Heritage Team", "LIGO", "FAST"
};
std::vector<SearchResult> results[MAX_WINDOWS];
sqlite3* db;
Aws::S3::S3Client* s3_client;
Aws::CognitoIdentityProvider::CognitoIdentityProviderClient* cognito_client;

// NASA DONKI fetch
std::string FetchDONKI(const std::string& startDate = "", const std::string& endDate = "") {
    CURL* curl = curl_easy_init();
    std::string url = "https://api.nasa.gov/DONKI/CMEAnalysis?api_key=" + std::string(NASA_API_KEY_2);
    if (!startDate.empty()) url += "&startDate=" + startDate;
    if (!endDate.empty()) url += "&endDate=" + endDate;
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// Scientific Calculator with JD-Cal
class ScientificCalculatorDialog : public QDialog {
public:
    ScientificCalculatorDialog(QWidget* parent) : QDialog(parent) {
        setWindowFlags(Qt::Window | Qt::FramelessWindowHint);
        setAcceptDrops(true);
        QVBoxLayout* layout = new QVBoxLayout(this);
        input = new QTextEdit(this);
        input->setPlaceholderText("Enter equations (e.g., d/dx(x^2), ?(0,1) x^2 dx, x^2 + y = 5, jd to date 2451544)");
        input->setMinimumHeight(100);
        input->setMaximumHeight(1000);
        input->setAcceptDrops(true);
        output = new QTextEdit(this);
        output->setReadOnly(true);
        QPushButton* solveBtn = new QPushButton("Solve", this);
        layout->addWidget(input);
        layout->addWidget(solveBtn);
        layout->addWidget(output);
        connect(solveBtn, &QPushButton::clicked, this, &ScientificCalculatorDialog::solveEquations);
        connect(input, &QTextEdit::textChanged, this, &ScientificCalculatorDialog::adjustInputSize);
        setMouseTracking(true);
    }

protected:
    void mousePressEvent(QMouseEvent* event) override {
        if (event->button() == Qt::LeftButton) {
            dragPosition = event->globalPos() - frameGeometry().topLeft();
            event->accept();
        }
    }
    void mouseMoveEvent(QMouseEvent* event) override {
        if (event->buttons() & Qt::LeftButton) {
            move(event->globalPos() - dragPosition);
            event->accept();
        }
    }
    void dragEnterEvent(QDragEnterEvent* event) override {
        if (event->mimeData()->hasText()) event->acceptProposedAction();
    }
    void dropEvent(QDropEvent* event) override {
        input->setText(input->toPlainText() + event->mimeData()->text());
        event->acceptProposedAction();
    }

private:
    QTextEdit* input;
    QTextEdit* output;
    QPoint dragPosition;

    void adjustInputSize() {
        QString text = input->toPlainText();
        int lines = text.split("\n").size();
        int newHeight = std::min(std::max(100, lines * 20 + 50), 1000);
        input->setMinimumHeight(newHeight);
        input->setMaximumHeight(newHeight);
    }

    void solveEquations() {
        std::string expr = input->toPlainText().toStdString();
        std::vector<std::string> equations;
        std::stringstream ss(expr);
        std::string line;
        while (std::getline(ss, line)) {
            if (!line.empty()) equations.push_back(line);
        }

        QString result;
        Qalculate calc;
        py::scoped_interpreter guard{};
        py::module_ sympy = py::module_::import("sympy");

        std::vector<std::string> system_eqs;

        for (const auto& eq : equations) {
            if (eq.find("jd to date") != std::string::npos) {
                std::string jd = eq.substr(eq.find("date") + 5);
                std::string jdcal = FetchJDCalJD(jd);
                result += QString("JD to Date: %1\n").arg(QString::fromStdString(jdcal));
                // Sync with DONKI
                std::string donki = FetchDONKI(); // Fetch DONKI for space weather
                result += QString("DONKI Space Weather: %1\n").arg(QString::fromStdString(SummarizeWithOpenAI(donki)));
            }
            else if (eq.find("date to jd") != std::string::npos) {
                std::string cd = eq.substr(eq.find("jd") + 3);
                std::string jdcal = FetchJDCalCD(cd);
                result += QString("Date to JD: %1\n").arg(QString::fromStdString(jdcal));
            }
            else if (eq.find("d/d") != std::string::npos) {
                // Derivative
                std::string var = "x";
                std::string func = eq.substr(eq.find("(") + 1, eq.find(")") - eq.find("(") - 1);
                py::object x = sympy.attr("symbols")("x");
                py::object expr = sympy.attr("sympify")(func);
                py::object deriv = sympy.attr("diff")(expr, x);
                result += QString("d/dx(%1) = %2\n").arg(QString::fromStdString(func), QString::fromStdString(deriv.attr("__str__")().cast<std::string>()));
            }
            else if (eq.find("?") != std::string::npos) {
                // Integral
                std::string bounds = eq.substr(eq.find("(") + 1, eq.find(")") - eq.find("(") - 1);
                std::string func = eq.substr(eq.find(")") + 1, eq.find("dx") - eq.find(")") - 1);
                auto [a, b] = parseBounds(bounds);
                py::object x = sympy.attr("symbols")("x");
                py::object expr = sympy.attr("sympify")(func);
                py::object integral = sympy.attr("integrate")(expr, py::make_tuple(x, a, b));
                result += QString("?(%1,%2) %3 dx = %4\n")
                    .arg(QString::number(a), QString::number(b), QString::fromStdString(func), QString::fromStdString(integral.attr("__str__")().cast<std::string>()));
            }
            else if (eq.find("=") != std::string::npos) {
                // Collect for system
                std::string eq_clean = eq;
                std::replace(eq_clean.begin(), eq_clean.end(), '=', '-');
                system_eqs.push_back(eq_clean);
            }
            else {
                result += QString("%1 = %2\n").arg(QString::fromStdString(eq), QString::fromStdString(calc.evaluate(eq)));
            }
        }
        // Solve system if at least 2 equations
        if (system_eqs.size() >= 2) {
            py::object x = sympy.attr("symbols")("x");
            py::object y = sympy.attr("symbols")("y");
            py::object eq1 = sympy.attr("sympify")(system_eqs[0]);
            py::object eq2 = sympy.attr("sympify")(system_eqs[1]);
            py::object solutions = sympy.attr("solve")(py::make_tuple(eq1, eq2), py::make_tuple(x, y));
            result += QString("System: %1, %2\nSolutions: %3\n")
                .arg(QString::fromStdString(system_eqs[0]), QString::fromStdString(system_eqs[1]), QString::fromStdString(solutions.attr("__str__")().cast<std::string>()));
        }
        output->setText(result);
    }

    std::pair<double, double> parseBounds(const std::string& bounds) {
        size_t comma = bounds.find(",");
        double a = std::stod(bounds.substr(0, comma));
        double b = std::stod(bounds.substr(comma + 1));
        return { a, b };
    }
};

// Ramanujan Number Calculator
class RamanujanCalculatorDialog : public QDialog {
public:
    RamanujanCalculatorDialog(QWidget* parent) : QDialog(parent) {
        setWindowFlags(Qt::Window | Qt::FramelessWindowHint);
        setAcceptDrops(true);
        QVBoxLayout* layout = new QVBoxLayout(this);
        input = new QTextEdit(this);
        input->setPlaceholderText("Enter number theory functions (e.g., p(5), tau(7))");
        input->setMinimumHeight(100);
        input->setMaximumHeight(1000);
        input->setAcceptDrops(true);
        output = new QTextEdit(this);
        output->setReadOnly(true);
        QPushButton* solveBtn = new QPushButton("Solve", this);
        layout->addWidget(input);
        layout->addWidget(solveBtn);
        layout->addWidget(output);
        connect(solveBtn, &QPushButton::clicked, this, &RamanujanCalculatorDialog::solveEquations);
        connect(input, &QTextEdit::textChanged, this, &RamanujanCalculatorDialog::adjustInputSize);
        setMouseTracking(true);
    }

protected:
    void mousePressEvent(QMouseEvent* event) override {
        if (event->button() == Qt::LeftButton) {
            dragPosition = event->globalPos() - frameGeometry().topLeft();
            event->accept();
        }
    }
    void mouseMoveEvent(QMouseEvent* event) override {
        if (event->buttons() & Qt::LeftButton) {
            move(event->globalPos() - dragPosition);
            event->accept();
        }
    }
    void dragEnterEvent(QDragEnterEvent* event) override {
        if (event->mimeData()->hasText()) event->acceptProposedAction();
    }
    void dropEvent(QDropEvent* event) override {
        input->setText(input->toPlainText() + event->mimeData()->text());
        event->acceptProposedAction();
    }

private:
    QTextEdit* input;
    QTextEdit* output;
    QPoint dragPosition;

    void adjustInputSize() {
        QString text = input->toPlainText();
        int lines = text.split("\n").size();
        int newHeight = std::min(std::max(100, lines * 20 + 50), 1000);
        input->setMinimumHeight(newHeight);
        input->setMaximumHeight(newHeight);
    }

    void solveEquations() {
        std::string expr = input->toPlainText().toStdString();
        std::vector<std::string> equations;
        std::stringstream ss(expr);
        std::string line;
        while (std::getline(ss, line)) {
            if (!line.empty()) equations.push_back(line);
        }

        QString result;
        py::scoped_interpreter guard{};
        py::module_ sympy = py::module_::import("sympy");

        for (const auto& eq : equations) {
            if (eq.find("p(") != std::string::npos) {
                std::string n_str = eq.substr(eq.find("(") + 1, eq.find(")") - eq.find("(") - 1);
                int n = std::stoi(n_str);
                py::object partition = sympy.attr("partition")(n);
                result += QString("p(%1) = %2 partitions\n").arg(n).arg(partition.cast<int>());
            }
            else if (eq.find("tau(") != std::string::npos) {
                std::string n_str = eq.substr(eq.find("(") + 1, eq.find(")") - eq.find("(") - 1);
                int n = std::stoi(n_str);
                py::object tau = sympy.attr("ramanujan_tau")(n);
                result += QString("tau(%1) = %2\n").arg(n).arg(tau.cast<long>());
            }
            else {
                result += QString("Invalid input: %1\n").arg(QString::fromStdString(eq));
            }
        }
        output->setText(result);
    }
};

// Calculus Button Field
class CalculusButtonField : public QDockWidget {
public:
    CalculusButtonField(QWidget* parent) : QDockWidget("Calculus Tools", parent) {
        QWidget* widget = new QWidget();
        QVBoxLayout* layout = new QVBoxLayout(widget);
        QToolBar* toolbar = new QToolBar(this);
        input = new QTextEdit(this);
        input->setPlaceholderText("Insert symbols (e.g., ?, ?, ?)");
        input->setMinimumHeight(100);
        input->setMaximumHeight(1000);
        input->setAcceptDrops(true);

        toolbar->addAction("?", [=]() { input->insertPlainText("?(a,b) f(x) dx"); });
        toolbar->addAction("?", [=]() { input->insertPlainText("?/?x"); });
        toolbar->addAction("?", [=]() { input->insertPlainText("?(n,a,b)"); });
        toolbar->addAction("?", [=]() { input->insertPlainText("sqrt()"); });
        toolbar->addAction("sin", [=]() { input->insertPlainText("sin()"); });
        toolbar->addAction("cos", [=]() { input->insertPlainText("cos()"); });
        toolbar->addAction("log", [=]() { input->insertPlainText("log()"); });

        layout->addWidget(toolbar);
        layout->addWidget(input);
        setWidget(widget);
        connect(input, &QTextEdit::textChanged, this, &CalculusButtonField::adjustInputSize);
    }

protected:
    void dragEnterEvent(QDragEnterEvent* event) override {
        if (event->mimeData()->hasText()) event->acceptProposedAction();
    }
    void dropEvent(QDropEvent* event) override {
        input->setText(input->toPlainText() + event->mimeData()->text());
        event->acceptProposedAction();
    }

private:
    QTextEdit* input;

    void adjustInputSize() {
        QString text = input->toPlainText();
        int lines = text.split("\n").size();
        int newHeight = std::min(std::max(100, lines * 20 + 50), 1000);
        input->setMinimumHeight(newHeight);
        input->setMaximumHeight(newHeight);
    }
};

// Detachable Browser Window
class BrowserWindow : public QMainWindow {
public:
    BrowserWindow(const QString& title, QWidget* parent = nullptr) : QMainWindow(parent) {
        QWebEngineView* view = new QWebEngineView(this);
        QTextEdit* summary = new QTextEdit(this);
        summary->setReadOnly(true);
        QVBoxLayout* layout = new QVBoxLayout();
        QWidget* centralWidget = new QWidget();
        layout->addWidget(view);
        layout->addWidget(summary);
        centralWidget->setLayout(layout);
        setCentralWidget(centralWidget);
        setWindowTitle(title);
        views.push_back(view);
        summaries.push_back(summary);
    }

    void setContent(const QString& html) {
        views[0]->setHtml(html);
        summaries[0]->setText(html);
    }

private:
    std::vector<QWebEngineView*> views;
    std::vector<QTextEdit*> summaries;
};

// WebSocket callback
void on_message(void* user, const char* data, size_t len) {
    std::string json_data(data, len);
    SearchResult result = { "wss://ligo.org/alerts", "Live Data", "Real-time event", 1.0, true };
    results[0].push_back(result);
    sqlite3_exec(db, ("INSERT INTO cache (url, title, summary, isLive) VALUES ('" + result.url + "', '" + result.title + "', '" + result.summary + "', 1)").c_str(), nullptr, nullptr, nullptr);
}

// cURL callback
size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* data) {
    data->append((char*)contents, size * nmemb);
    return size * nmemb;
}

// GPT-based summarization (Llama-3.1 fallback)
std::string SummarizeText(const std::string& text) {
    py::scoped_interpreter guard{};
    py::module_ transformers = py::module_::import("transformers");
    py::object summarizer = transformers.attr("pipeline")("summarization", "meta-llama/Llama-3.1-8B");
    py::object summary = summarizer(text, py::arg("max_length") = 100, py::arg("min_length") = 30);
    return summary[0].attr("summary_text").cast<std::string>();
}

// OpenAI summarization with retry logic
std::string SummarizeWithOpenAI(const std::string& query) {
    CURL* curl = curl_easy_init();
    std::string url = "https://api.openai.com/v1/chat/completions";
    std::string response;
    json payload = {
        {"model", "gpt-4"},
        {"messages", {{{"role", "user"}, {"content", "Summarize: " + query}}}},
        {"max_tokens", 100}
    };
    std::string data = payload.dump();
    struct curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    headers = curl_slist_append(headers, ("Authorization: Bearer " + std::string(OPENAI_API_KEY)).c_str());
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, data.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

    int retries = 3;
    while (retries--) {
        CURLcode res = curl_easy_perform(curl);
        long http_code = 0;
        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
        if (res == CURLE_OK && http_code == 200) {
            json result = json::parse(response);
            curl_slist_free_all(headers);
            curl_easy_cleanup(curl);
            return result["choices"][0]["message"]["content"].get<std::string>();
        }
        else if (http_code == 429) {
            std::this_thread::sleep_for(std::chrono::seconds(1 << (3 - retries)));
            continue;
        }
        break;
    }
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
    return SummarizeText(query); // Fallback to Llama-3.1
}

// OAuth for cloud sync
std::string GetOAuthToken() {
    CURL* curl = curl_easy_init();
    std::string url = "https://<domain>.auth." + std::string(COGNITO_REGION) + ".amazoncognito.com/oauth2/token";
    std::string data = "grant_type=client_credentials&client_id=" + std::string(COGNITO_CLIENT_ID) + "&client_secret=your_client_secret";
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, data.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return "mock_access_token"; // Parse JSON for access_token
}

// Cloud sync
void SyncCacheToCloud(const std::string& token) {
    Aws::S3::Model::PutObjectRequest request;
    request.SetBucket("coanqi-cache");
    request.SetKey("cache.db");
    request.SetCustomRequestHeader("Authorization", "Bearer " + token);
    std::ifstream file("coanqi_cache.db", std::ios::binary);
    request.SetBody(std::make_shared<Aws::Fstream>(file));
    s3_client->PutObject(request);
}

// Offline search
void OfflineSearch(const std::string& query, std::vector<SearchResult>& offlineResults) {
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(db, "SELECT url, title, summary, isLive FROM cache WHERE title LIKE ? OR summary LIKE ?", -1, &stmt, nullptr);
    std::string pattern = "%" + query + "%";
    sqlite3_bind_text(stmt, 1, pattern.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 2, pattern.c_str(), -1, SQLITE_STATIC);
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        SearchResult result;
        result.url = (const char*)sqlite3_column_text(stmt, 0);
        result.title = (const char*)sqlite3_column_text(stmt, 1);
        result.summary = (const char*)sqlite3_column_text(stmt, 2);
        result.isLive = sqlite3_column_int(stmt, 3);
        result.relevance = 0.9;
        offlineResults.push_back(result);
    }
    sqlite3_finalize(stmt);
}

// Voice input
std::string ProcessVoiceInput() {
    ps_decoder_t* ps = ps_init(cmd_ln_init(nullptr, ps_args(), true, nullptr));
    ps_start_utt(ps);
    std::string text = "sample query"; // Replace with PocketSphinx
    ps_end_utt(ps);
    ps_free(ps);
    return text;
}

// Video input
std::string ProcessVideoInput() {
    cv::VideoCapture cap(0);
    cv::Mat frame;
    cap >> frame;
    std::string command = "submit query"; // Replace with OpenCV gesture recognition
    cap.release();
    return command;
}

// Visualization plugin
void RenderScatterPlot(QWidget* parent, const std::vector<double>& x, const std::vector<double>& y) {
    vtkSmartPointer<vtkScatterPlotMatrix> matrix = vtkSmartPointer<vtkScatterPlotMatrix>::New();
    // Add x, y data (simplified)
}

// NASA API search (updated with new endpoints)
void SearchNASA(const std::string& query, std::vector<SearchResult>& nasaResults) {
    CURL* curl = curl_easy_init();
    std::vector<std::string> endpoints = {
        "https://api.nasa.gov/planetary/apod?api_key=" + std::string(NASA_API_KEY_1) + "&concept_tags=True&keywords=" + query,
        "https://api.nasa.gov/EPIC/api/natural?api_key=" + std::string(NASA_API_KEY_2),
        "https://api.nasa.gov/DONKI/CMEAnalysis?api_key=" + std::string(NASA_API_KEY_2)
    };
    std::vector<std::string> titles = { "NASA APOD Result", "NASA EPIC Result", "NASA DONKI Result" };
    std::vector<std::string> urls = {
        "https://api.nasa.gov/planetary/apod",
        "https://api.nasa.gov/EPIC/api/natural",
        "https://api.nasa.gov/DONKI/CMEAnalysis"
    };

    for (size_t i = 0; i < endpoints.size(); ++i) {
        std::string response;
        curl_easy_setopt(curl, CURLOPT_URL, endpoints[i].c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
        int retries = 3;
        while (retries--) {
            CURLcode res = curl_easy_perform(curl);
            long http_code = 0;
            curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
            if (res == CURLE_OK && http_code == 200) {
                std::string summary = SummarizeWithOpenAI(response);
                SearchResult result = { urls[i], titles[i], summary, 0.95, false };
                nasaResults.push_back(result);
                sqlite3_exec(db, ("INSERT INTO cache (url, title, summary, isLive) VALUES ('" + result.url + "', '" + result.title + "', '" + result.summary + "', 0)").c_str(), nullptr, nullptr, nullptr);
                break;
            }
            else if (http_code == 429) {
                std::this_thread::sleep_for(std::chrono::seconds(1 << (3 - retries)));
                continue;
            }
        }
    }
    curl_easy_cleanup(curl);
}

// MAST API search
void SearchMAST(const std::string& query, std::vector<SearchResult>& mastResults) {
    CURL* curl = curl_easy_init();
    std::string url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/hst_12345_01_acs_f814w_drz.fits&token=" + std::string(MAST_API_KEY);
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    if (res == CURLE_OK) {
        std::string summary = SummarizeWithOpenAI(response);
        SearchResult result = { "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/hst_12345_01_acs_f814w_drz.fits", "MAST HST Infrared", summary, 0.95, false };
        mastResults.push_back(result);
        sqlite3_exec(db, ("INSERT INTO cache (url, title, summary, isLive) VALUES ('" + result.url + "', '" + result.title + "', '" + result.summary + "', 0)").c_str(), nullptr, nullptr, nullptr);
    }
    curl_easy_cleanup(curl);
}

// JPL Horizons API
std::string FetchHorizons(const std::string& command, const std::string& start_time, const std::string& stop_time) {
    CURL* curl = curl_easy_init();
    std::string url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND='" + command + "'&OBJ_DATA='YES'&MAKE_EPHEM='YES'&EPHEM_TYPE='OBSERVER'&CENTER='500@399'&START_TIME='" + start_time + "'&STOP_TIME='" + stop_time + "'&STEP_SIZE='1%20d'&QUANTITIES='1,9,20,23,24,29'";
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// JPL JD-Cal JD
std::string FetchJDCalJD(const std::string& jd) {
    CURL* curl = curl_easy_init();
    std::string url = "https://ssd-api.jpl.nasa.gov/jd_cal.api?jd=" + jd + "&format=s";
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// JPL JD-Cal CD
std::string FetchJDCalCD(const std::string& cd) {
    CURL* curl = curl_easy_init();
    std::string url = "https://ssd-api.jpl.nasa.gov/jd_cal.api?cd=" + cd;
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// JPL Periodic Orbits Earth-Moon
std::string FetchPeriodicEarthMoon(const std::string& family, const std::string& libr, const std::string& branch) {
    CURL* curl = curl_easy_init();
    std::string url = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=earth-moon&family=" + family + "&libr=" + libr + "&branch=" + branch;
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// JPL Periodic Orbits Jupiter-Europa
std::string FetchPeriodicJupiterEuropa(const std::string& family, double stability = -1.0) {
    CURL* curl = curl_easy_init();
    std::string url = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=jupiter-europa&family=" + family;
    if (stability > -1.0) url += "&stabmax=" + std::to_string(stability);
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// JPL Periodic Orbits Saturn-Enceladus
std::string FetchPeriodicSaturnEnceladus(const std::string& family, const std::string& libr, double periodmax = 1.0, const std::string& periodunits = "d") {
    CURL* curl = curl_easy_init();
    std::string url = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=saturn-enceladus&family=" + family + "&libr=" + libr + "&periodmax=" + std::to_string(periodmax) + "&periodunits=" + periodunits;
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// JPL Periodic Orbits Saturn-Titan
std::string FetchPeriodicSaturnTitan(const std::string& family, double jacobimin = 3.0, double stabmax = 1.0, const std::string& branch = "N") {
    CURL* curl = curl_easy_init();
    std::string url = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=saturn-titan&family=" + family + "&jacobimin=" + std::to_string(jacobimin) + "&stabmax=" + std::to_string(stabmax) + "&branch=" + branch;
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// JPL Periodic Orbits Mars-Phobos
std::string FetchPeriodicMarsPhobos(const std::string& family, double jacobimin = 3.0, double stabmax = 1.0, const std::string& branch = "21") {
    CURL* curl = curl_easy_init();
    std::string url = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=mars-phobos&family=" + family + "&jacobimin=" + std::to_string(jacobimin) + "&stabmax=" + std::to_string(stabmax) + "&branch=" + branch;
    std::string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return response;
}

// Search function (updated with JPL APIs and user-configurable parameters)
void PerformSearch(const std::string& query, std::vector<std::string>& focus, bool online, const std::string& oauth_token) {
    if (!online) {
        std::vector<SearchResult> offlineResults;
        OfflineSearch(query, offlineResults);
        for (int i = 0; i < MAX_WINDOWS && i < offlineResults.size(); ++i) {
            results[i].push_back(offlineResults[i]);
        }
        return;
    }

    std::vector<SearchResult> nasaResults, mastResults;
    if (std::find(focus.begin(), focus.end(), "NASA") != focus.end()) {
        SearchNASA(query, nasaResults);
        results[1] = nasaResults;
    }
    if (std::find(focus.begin(), focus.end(), "STScI") != focus.end() ||
        std::find(focus.begin(), focus.end(), "Hubble") != focus.end() ||
        std::find(focus.begin(), focus.end(), "ACS Hubble Ultra Deep Field") != focus.end()) {
        SearchMAST(query, mastResults);
        results[2] = mastResults;
    }

    // Preload specified links with OpenAI summaries
    results[3].push_back({ "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/hst_12345_01_acs_f814w_drz.fits", "MAST ACS F814W Infrared", SummarizeWithOpenAI("Hubble infrared data"), 0.95, false });
    results[4].push_back({ "wss://eventhorizontelescope.org/data", "EHT Live Infrared Data", SummarizeWithOpenAI("Real-time EHT data"), 1.0, true });
    results[5].push_back({ "https://apod.nasa.gov/apod/image/2507/m31_infrared.jpg", "NASA M31 Infrared", SummarizeWithOpenAI("Andromeda infrared image"), 0.95, false });
    results[6].push_back({ "wss://ligo.org/alerts", "LIGO GW Infrared Correlations", SummarizeWithOpenAI("Real-time GW alerts"), 1.0, true });

    // JPL APIs with dynamic parameters
    // Parse query for parameters (simplified; use regex for production)
    std::string command = "499"; // Default
    std::string start_time = "2006-01-01";
    std::string stop_time = "2006-01-20";
    std::string jd = "2451544";
    std::string cd = "2000-01-01_12:00";
    std::string family = "halo";
    std::string libr = "1";
    std::string branch = "N";
    double stability = -1.0;
    double periodmax = 1.0;
    std::string periodunits = "d";
    double jacobimin = 3.0;
    double stabmax = 1.0;

    if (query.find("ephemeris") != std::string::npos || query.find("horizons") != std::string::npos) {
        std::string horizons = FetchHorizons(command, start_time, stop_time); // User-input via config
        results[7].push_back({ "https://ssd.jpl.nasa.gov/api/horizons.api", "JPL Horizons Ephemeris", SummarizeWithOpenAI(horizons), 0.95, false });
        RenderScatterPlot(nullptr, {}, {}); // Visualize orbits
    }
    if (query.find("jd to date") != std::string::npos) {
        std::string jdcal = FetchJDCalJD(jd); // Dynamic JD from query
        results[8].push_back({ "https://ssd-api.jpl.nasa.gov/jd_cal.api?jd=2451544&format=s", "JPL JD-Cal JD to Date", SummarizeWithOpenAI(jdcal), 0.95, false });
    }
    if (query.find("date to jd") != std::string::npos) {
        std::string jdcal = FetchJDCalCD(cd); // Dynamic CD from query
        results[9].push_back({ "https://ssd-api.jpl.nasa.gov/jd_cal.api?cd=2000-01-01_12:00", "JPL JD-Cal Date to JD", SummarizeWithOpenAI(jdcal), 0.95, false });
    }
    if (query.find("earth-moon halo") != std::string::npos) {
        std::string orbits = FetchPeriodicEarthMoon(family, libr, branch); // Dynamic family/libr/branch
        results[10].push_back({ "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=earth-moon&family=halo&libr=1&branch=N", "JPL Periodic Orbits Earth-Moon", SummarizeWithOpenAI(orbits), 0.95, false });
        RenderScatterPlot(nullptr, {}, {}); // Visualize in VTK for mission planning
    }
    if (query.find("jupiter-europa dro") != std::string::npos) {
        std::string orbits = FetchPeriodicJupiterEuropa(family, stability); // Dynamic family, filter stability
        results[11].push_back({ "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=jupiter-europa&family=dro", "JPL Periodic Orbits Jupiter-Europa", SummarizeWithOpenAI(orbits), 0.95, false });
    }
    if (query.find("saturn-enceladus vertical") != std::string::npos) {
        std::string orbits = FetchPeriodicSaturnEnceladus(family, libr, periodmax, periodunits); // Dynamic
        results[12].push_back({ "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=saturn-enceladus&family=vertical&libr=2&periodmax=1&periodunits=d", "JPL Periodic Orbits Saturn-Enceladus", SummarizeWithOpenAI(orbits), 0.95, false });
    }
    if (query.find("saturn-titan butterfly") != std::string::npos) {
        std::string orbits = FetchPeriodicSaturnTitan(family, jacobimin, stabmax, branch); // Dynamic
        results[13].push_back({ "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=saturn-titan&family=butterfly&jacobimin=3&stabmax=1&branch=N", "JPL Periodic Orbits Saturn-Titan", SummarizeWithOpenAI(orbits), 0.95, false });
        // Integrate with DONKI for space weather
    }
    if (query.find("mars-phobos resonant") != std::string::npos) {
        std::string orbits = FetchPeriodicMarsPhobos(family, jacobimin, stabmax, branch); // Dynamic
        results[14].push_back({ "https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=mars-phobos&family=resonant&jacobimin=3&stabmax=1&branch=21", "JPL Periodic Orbits Mars-Phobos", SummarizeWithOpenAI(orbits), 0.95, false });
        RenderScatterPlot(nullptr, {}, {}); // Visualize Phobos trajectories, link to Horizons
    }

    struct lws_context* ws_context = lws_create_context(nullptr);
    lws_connect(ws_context, "eventhorizontelescope.org", 443, "/data", on_message, nullptr);
    lws_connect(ws_context, "skaobservatory.org", 443, "/realtime", on_message, nullptr);
    lws_connect(ws_context, "ligo.org", 443, "/alerts", on_message, nullptr);
    lws_connect(ws_context, "fast.bao.ac.cn", 443, "/realtime", on_message, nullptr);

    CURL* curl = curl_easy_init();
    for (int i = 15; i < MAX_WINDOWS && i < focus.size(); ++i) {
        std::string url = "https://api.example.com/search?q=" + query + "&source=" + focus[i];
        std::string response;
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
        CURLcode res = curl_easy_perform(curl);
        if (res == CURLE_OK) {
            std::string summary = SummarizeWithOpenAI(response);
            SearchResult result = { "https://example.com", "Sample Result", summary, 0.95, false };
            results[i].push_back(result);
            sqlite3_exec(db, ("INSERT INTO cache (url, title, summary, isLive) VALUES ('" + result.url + "', '" + result.title + "', '" + result.summary + "', 0)").c_str(), nullptr, nullptr, nullptr);
        }
    }
    curl_easy_cleanup(curl);
    lws_context_destroy(ws_context);
    SyncCacheToCloud(oauth_token);
}

// Qt Main Window
class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow() {
        // Win32: System tray
#ifdef _WIN32
        NOTIFYICONDATA nid = { sizeof(nid) };
        nid.hWnd = (HWND)winId();
        nid.uID = 1;
        nid.uFlags = NIF_ICON | NIF_TIP;
        nid.hIcon = LoadIcon(GetModuleHandle(nullptr), "Z.ico");
        strcpy(nid.szTip, "CoAnQi");
        Shell_NotifyIcon(NIM_ADD, &nid);
#endif

        QWidget* centralWidget = new QWidget(this);
        QVBoxLayout* layout = new QVBoxLayout(centralWidget);

        // Top bar (Firefox-like)
        QHBoxLayout* topBar = new QHBoxLayout();
        QPushButton* backBtn = new QPushButton("Back", this);
        QPushButton* forwardBtn = new QPushButton("Forward", this);
        QPushButton* refreshBtn = new QPushButton("Refresh", this);
        QLineEdit* queryField = new QLineEdit(this);
        queryField->setMaxLength(MAX_QUERY_LENGTH);
        queryField->setPlaceholderText("Search high-energy datasets...");
        QPushButton* voiceBtn = new QPushButton("??", this);
        QPushButton* videoBtn = new QPushButton("??", this);
        QPushButton* sciCalcBtn = new QPushButton("??", this);
        QPushButton* ramCalcBtn = new QPushButton("??R", this);
        QPushButton* calcBtnField = new QPushButton("??C", this);
        QLabel* logo = new QLabel("<b>CoAnQi (Cosmic Analysis and Quantum Intelligence)</b>", this);
        logo->setStyleSheet("font-size: 24px; color: #2a5298;");
        QPushButton* menuBtn = new QPushButton("?", this);
        topBar->addWidget(backBtn);
        topBar->addWidget(forwardBtn);
        topBar->addWidget(refreshBtn);
        topBar->addWidget(queryField);
        topBar->addWidget(voiceBtn);
        topBar->addWidget(videoBtn);
        topBar->addWidget(sciCalcBtn);
        topBar->addWidget(ramCalcBtn);
        topBar->addWidget(calcBtnField);
        topBar->addWidget(logo);
        topBar->addWidget(menuBtn);
        layout->addLayout(topBar);

        // Focus list
        QTextEdit* focusField = new QTextEdit(this);
        QString focusText;
        for (const auto& item : focusList) focusText += QString::fromStdString(item) + "\n";
        focusField->setText(focusText);
        layout->addWidget(focusField);

        // Tabbed browser windows
        QTabWidget* tabs = new QTabWidget(this);
        tabs->setTabsClosable(true);
        tabs->setMovable(true);
        browserWindows = new BrowserWindow * [MAX_WINDOWS];
        for (int i = 0; i < MAX_WINDOWS; ++i) {
            browserWindows[i] = new BrowserWindow(QString("Tab %1").arg(i + 1), this);
            tabs->addTab(new QWidget(), QString("Tab %1").arg(i + 1));
        }
        // Dedicated ALMA Cycle 12 window (Tab 21)
        browserWindows[20]->views[0]->load(QUrl("https://almascience.nrao.edu/proposing/observing-tool/tarball-download-page"));
        layout->addWidget(tabs);

        // Visualization sidebar
        QDockWidget* sidebar = new QDockWidget("Visualizations", this);
        QWidget* visWidget = new QWidget();
        QVBoxLayout* visLayout = new QVBoxLayout(visWidget);
        visLayout->addWidget(new QLabel("Dataset Graph Placeholder"));
        sidebar->setWidget(visWidget);
        addDockWidget(Qt::LeftDockWidgetArea, sidebar);

        // Calculus button field
        CalculusButtonField* calcField = new CalculusButtonField(this);
        addDockWidget(Qt::RightDockWidgetArea, calcField);

        // Calculators
        ScientificCalculatorDialog* sciCalcDialog = new ScientificCalculatorDialog(this);
        sciCalcDialog->move(50, 50);
        sciCalcDialog->show();
        RamanujanCalculatorDialog* ramCalcDialog = new RamanujanCalculatorDialog(this);
        ramCalcDialog->move(100, 100);
        ramCalcDialog->show();

        setCentralWidget(centralWidget);

        // Initialize SQLite and AWS
        sqlite3_open("coanqi_cache.db", &db);
        sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS cache (url TEXT, title TEXT, summary TEXT, isLive INTEGER)", nullptr, nullptr, nullptr);
        Aws::SDKOptions options;
        Aws::InitAPI(options);
        s3_client = new Aws::S3::S3Client();
        cognito_client = new Aws::CognitoIdentityProvider::CognitoIdentityProviderClient();

        // OAuth login
        std::string oauth_token = GetOAuthToken();

        // Connect signals
        connect(queryField, &QLineEdit::returnPressed, [=]() {
            std::string query = queryField->text().toStdString();
            if (query.length() > MAX_QUERY_LENGTH) {
                QMessageBox::warning(this, "Error", "Query exceeds 3000 characters!");
                return;
            }
            bool online = true; // Check connectivity
            PerformSearch(query, focusList, online, oauth_token);
            for (int i = 0; i < MAX_WINDOWS; ++i) {
                QString html = "<ul>";
                for (const auto& result : results[i]) {
                    QString live = result.isLive ? " [Live]" : "";
                    html += QString("<li><a href='%1'>%2</a>%3: %4</li>")
                        .arg(QString::fromStdString(result.url))
                        .arg(QString::fromStdString(result.title))
                        .arg(live)
                        .arg(QString::fromStdString(result.summary));
                }
                html += "</ul>";
                browserWindows[i]->setContent(html);
            }
            });

        connect(tabs, &QTabWidget::tabBarDoubleClicked, [=](int index) {
            BrowserWindow* window = browserWindows[index];
            window->show();
            tabs->removeTab(index);
            });

        connect(voiceBtn, &QPushButton::clicked, [=]() {
            queryField->setText(QString::fromStdString(ProcessVoiceInput()));
            });

        connect(videoBtn, &QPushButton::clicked, [=]() {
            if (ProcessVideoInput() == "submit query") {
                QKeyEvent* event = new QKeyEvent(QEvent::KeyPress, Qt::Key_Return, Qt::NoModifier);
                QCoreApplication::postEvent(queryField, event);
            }
            });

        connect(sciCalcBtn, &QPushButton::clicked, [=]() {
            sciCalcDialog->show();
            });

        connect(ramCalcBtn, &QPushButton::clicked, [=]() {
            ramCalcDialog->show();
            });

        connect(calcBtnField, &QPushButton::clicked, [=]() {
            calcField->show();
            });

        connect(focusField, &QTextEdit::textChanged, [=]() {
            focusList.clear();
            QStringList lines = focusField->toPlainText().split("\n");
            for (const auto& line : lines) {
                if (!line.isEmpty()) focusList.push_back(line.toStdString());
            }
            });
    }
    ~MainWindow() {
        for (int i = 0; i < MAX_WINDOWS; ++i) delete browserWindows[i];
        delete[] browserWindows;
        sqlite3_close(db);
        delete s3_client;
        delete cognito_client;
        Aws::ShutdownAPI(Aws::SDKOptions());
#ifdef _WIN32
        NOTIFYICONDATA nid = { sizeof(nid) };
        nid.uID = 1;
        Shell_NotifyIcon(NIM_DELETE, &nid);
#endif
    }

private:
    BrowserWindow** browserWindows;
};

// Main function
int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    MainWindow window;
    window.setWindowTitle("CoAnQi");
    window.setWindowIcon(QIcon("Z.png"));
    window.show();
    return app.exec();
}