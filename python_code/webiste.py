import sys
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QGuiApplication, QFont
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QTextEdit, QLabel, QHBoxLayout


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Peak noPeak identification Website")
        self.setGeometry(100, 100, 800, 600)

        self.total_pdf_files = 100
        self.files_per_page = 10
        self.current_page = 1
        self.total_pages = (self.total_pdf_files + self.files_per_page - 1) // self.files_per_page

        # Create buttons for PDF files
        self.pdf_buttons_layout = QVBoxLayout()
        pdf_buttons_widget = QWidget()
        pdf_buttons_widget.setLayout(self.pdf_buttons_layout)

        self.pdf_buttons = []
        self.create_pdf_buttons()

        self.setCentralWidget(pdf_buttons_widget)

        # Create navigation buttons
        navigation_layout = QHBoxLayout()
        navigation_widget = QWidget()
        navigation_widget.setLayout(navigation_layout)

        prev_button = QPushButton("Previous")
        prev_button.clicked.connect(self.prev_page)
        navigation_layout.addWidget(prev_button)

        next_button = QPushButton("Next")
        next_button.clicked.connect(self.next_page)
        navigation_layout.addWidget(next_button)

        self.pdf_buttons_layout.addWidget(navigation_widget)

        # Initialize graph window
        self.graph_window = GraphWindow()

    def create_pdf_buttons(self):
        start_index = (self.current_page - 1) * self.files_per_page + 1
        end_index = min(start_index + self.files_per_page, self.total_pdf_files + 1)

        for i in range(start_index, end_index):
            pdf_button = QPushButton(f"PDF {i}")
            pdf_button.setStyleSheet("QPushButton { font-size: 18px; padding: 10px; }")
            pdf_button.clicked.connect(lambda _, idx=i: self.open_pdf(idx))
            self.pdf_buttons.append(pdf_button)
            self.pdf_buttons_layout.addWidget(pdf_button)

    def update_pdf_buttons(self):
        for button in self.pdf_buttons:
            self.pdf_buttons_layout.removeWidget(button)
            button.deleteLater()
        self.pdf_buttons = []
        self.create_pdf_buttons()

    def open_pdf(self, idx):
        self.graph_window.set_pdf_index(idx)
        self.graph_window.show()

    def prev_page(self):
        if self.current_page > 1:
            self.current_page -= 1
            self.update_pdf_buttons()

    def next_page(self):
        if self.current_page < self.total_pages:
            self.current_page += 1
            self.update_pdf_buttons()


class GraphWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Graph Window")
        self.setGeometry(200, 200, 800, 600)

        # Create buttons Y and N
        y_button = QPushButton("Y", self)
        y_button.setGeometry(20, 20, 120, 50)
        y_button.setStyleSheet("QPushButton { font-size: 22px; padding: 10px; background-color: green; color: white; }")
        y_button.clicked.connect(self.update_table)

        n_button = QPushButton("N", self)
        n_button.setGeometry(160, 20, 120, 50)
        n_button.setStyleSheet("QPushButton { font-size: 22px; padding: 10px; background-color: red; color: white; }")
        n_button.clicked.connect(self.revert_table)

        # Create graph display
        graph_display = QPushButton("Graph", self)
        graph_display.setGeometry(20, 100, 120, 50)
        graph_display.setStyleSheet("QPushButton { font-size: 22px; padding: 10px; background-color: #007bff; color: white; }")
        graph_display.clicked.connect(self.open_3d_graph)

        # Create table display for content
        self.table_display = QTextEdit(self)
        self.table_display.setGeometry(20, 180, 760, 380)
        self.table_display.setStyleSheet("QTextEdit { font-size: 18px; }")

    def set_pdf_index(self, idx):
        self.pdf_index = idx
        self.setWindowTitle(f"Graph Window - PDF {idx}")

    def update_table(self):
        content = f"<h2>PDF {self.pdf_index} - Y Button Clicked</h2>"
        self.table_display.setHtml(content)

    def revert_table(self):
        content = f"<h2>PDF {self.pdf_index} - N Button Clicked</h2>"
        self.table_display.setHtml(content)

    def open_3d_graph(self):
        # Code to open the graph in 3D
        pass


def main():
    app = QApplication(sys.argv)
    app.setFont(QFont("Arial", 12))

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()
