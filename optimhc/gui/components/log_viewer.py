"""
Log viewer component for optiMHC GUI.
"""

import os
from typing import List, Optional
import streamlit as st


def display_logs(logs: List[str]):
    """
    Display logs as read-only text with auto-scrolling.

    Args:
        logs: List of log messages
    """
    if not logs:
        st.info("No logs to display yet...")
        return

    # Join logs with newlines
    log_text = "\n".join(logs)

    # Use a container with custom CSS to create a taller scrollable area
    log_container = st.container()

    with log_container:
        # Add custom CSS for taller log area with scrollbar
        st.markdown("""
        <style>
        .stCodeBlock {
            max-height: 600px !important;
            overflow-y: auto !important;
            margin-bottom: 10px !important;
        }
        </style>
        """, unsafe_allow_html=True)

        # Use st.code for read-only display
        st.code(log_text, language="plain")

    # Add JavaScript to auto-scroll to bottom
    # This ensures the latest logs are always visible
    auto_scroll_js = """
    <script>
    // Method 1: Direct scrolling with timeout
    setTimeout(function() {
        const codeBlocks = parent.document.querySelectorAll('.stCodeBlock pre');
        const lastCodeBlock = codeBlocks[codeBlocks.length - 1];
        if (lastCodeBlock) {
            lastCodeBlock.scrollTop = lastCodeBlock.scrollHeight;
        }
    }, 100);

    // Method 2: MutationObserver to detect changes and scroll
    const observer = new MutationObserver((mutations) => {
        mutations.forEach((mutation) => {
            if (mutation.type === 'childList') {
                const codeBlocks = parent.document.querySelectorAll('.stCodeBlock pre');
                const lastCodeBlock = codeBlocks[codeBlocks.length - 1];
                if (lastCodeBlock) {
                    lastCodeBlock.scrollTop = lastCodeBlock.scrollHeight;
                }
            }
        });
    });

    // Observe code block for changes
    setTimeout(function() {
        const codeBlocks = parent.document.querySelectorAll('.stCodeBlock pre');
        const lastCodeBlock = codeBlocks[codeBlocks.length - 1];
        if (lastCodeBlock) {
            observer.observe(lastCodeBlock, { childList: true, subtree: true });
        }
    }, 200);
    </script>
    """
    st.components.v1.html(auto_scroll_js, height=0)


def find_pipeline_log_file() -> Optional[str]:
    """
    Find the pipeline log file based on configuration.

    Returns:
        Path to the log file or None if not found
    """
    # Try to get output directory and experiment name from config
    output_dir = None
    experiment_name = None

    if "config" in st.session_state:
        config = st.session_state.config
        output_dir = config.get("outputDir")
        experiment_name = config.get("experimentName")

    # First, check if we already know the log file path from pipeline execution
    if "pipeline_log_path" in st.session_state and st.session_state.pipeline_log_path:
        log_path = st.session_state.pipeline_log_path
        if os.path.exists(log_path):
            return log_path

    # Next, try to find log file in the expected pipeline output directory
    if output_dir and experiment_name:
        experiment_dir = os.path.join(output_dir, experiment_name)
        if os.path.exists(experiment_dir):
            # Pipeline's default log file
            log_path = os.path.join(experiment_dir, "log")
            if os.path.exists(log_path):
                return log_path

            # Look for any log files in the experiment directory
            for root, _, files in os.walk(experiment_dir):
                for file in files:
                    if file.endswith(".log") or file == "log":
                        return os.path.join(root, file)

    # If no log in experiment directory, search the main output directory
    if output_dir and os.path.exists(output_dir):
        log_files = []
        for root, _, files in os.walk(output_dir):
            log_files.extend([os.path.join(root, f) for f in files
                            if f.endswith(".log") or f == "log"])

        if log_files:
            # Return the most recently modified log file
            return max(log_files, key=os.path.getmtime)

    return None


def read_log_file(log_path, max_lines=1000):
    """
    Read log content from file.

    Args:
        log_path: Path to the log file
        max_lines: Maximum number of lines to read

    Returns:
        List of log lines
    """
    try:
        if not os.path.exists(log_path):
            return []

        with open(log_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            # Return the last max_lines lines
            return [line.rstrip() for line in lines[-max_lines:]]
    except Exception as e:
        print(f"Error reading log file: {str(e)}")
        return []


def update_logs():
    """
    Update logs from the log file.

    Returns:
        True if logs were updated, False otherwise
    """
    # Find the log file
    log_path = find_pipeline_log_file()
    if not log_path:
        return False

    # Read the log file
    logs = read_log_file(log_path)
    if not logs:
        return False

    # Update session state
    st.session_state.logs = logs
    return True


def log_viewer(process=None):
    """
    Simple log viewer with manual refresh button.

    Args:
        process: Optional subprocess to monitor for status
    """
    st.subheader("Log Output")

    # Status indicator if process is provided
    if process:
        if process.poll() is None:
            st.caption("📋 Process is running...")
        else:
            ret_code = process.poll()
            if ret_code == 0:
                st.caption("✅ Process completed successfully. Return code: 0")
            else:
                st.caption(f"⚠️ Process completed with errors. Return code: {ret_code}")

    # Find log file
    log_path = find_pipeline_log_file()

    # Controls row
    col1, col2, col3 = st.columns([2, 1, 1])

    with col1:
        if log_path:
            st.caption(f"Log file: {log_path}")
        else:
            st.caption("No log file found")

    with col2:
        # Clear logs button
        if st.button("Clear Logs"):
            st.session_state.logs = []
            st.rerun()

    with col3:
        # Refresh button
        if st.button("Refresh Logs"):
            update_logs()
            st.rerun()

    # Debug info (collapsed)
    with st.expander("Debug Info", expanded=False):
        log_path = find_pipeline_log_file() or "Not found"
        log_exists = "Yes" if log_path != "Not found" and os.path.exists(log_path) else "No"
        log_size = "0" if log_path == "Not found" or not os.path.exists(log_path) else str(os.path.getsize(log_path))

        process_info = ""
        if process:
            process_info = f"""
Process PID: {process.pid}
Process Return Code: {process.poll()}
Has stdout: {"Yes" if hasattr(process, 'stdout') and process.stdout else "No"}
"""

        st.code(f"""{process_info}
Log Count: {len(st.session_state.logs) if "logs" in st.session_state else 0}
Log File: {log_path}
Log File Exists: {log_exists}
Log File Size: {log_size} bytes
        """)

    # Initialize logs if needed
    if "logs" not in st.session_state:
        st.session_state.logs = []
        # Try to load logs the first time
        update_logs()

    # Display the logs
    display_logs(st.session_state.logs)
