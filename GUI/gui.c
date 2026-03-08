#include <gtk/gtk.h>
#include <signal.h>
#include <sys/wait.h>

typedef struct {
  const char *key;
  const char *label;
  GtkWidget *entry;
} ParamField;

typedef struct {
  GtkWidget *window;
  GtkWidget *entry_project;
  GtkWidget *entry_run_dir;
  GtkWidget *entry_bin;
  GtkWidget *spin_nproc;
  GtkWidget *combo_mode;
  GtkWidget *text_input;
  GtkWidget *text_log;
  GtkWidget *status_label;
  GtkWidget *btn_run;
  GtkWidget *btn_stop;
  GtkWidget *entry_plot_cmd;

  ParamField fields[19];
  int nfields;

  GPid child_pid;
  GIOChannel *stdout_ch;
  GIOChannel *stderr_ch;
  guint stdout_watch_id;
  guint stderr_watch_id;
  gboolean running;
} App;

typedef struct {
  GtkWidget *entry;
  GtkWidget *window;
  gboolean select_folder;
} PathPickerData;

static const int k_mode_values[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

static void append_log(App *app, const char *text) {
  GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(app->text_log));
  GtkTextIter end;
  gtk_text_buffer_get_end_iter(buf, &end);
  gtk_text_buffer_insert(buf, &end, text, -1);
  gtk_text_buffer_get_end_iter(buf, &end);
  gtk_text_view_scroll_to_iter(GTK_TEXT_VIEW(app->text_log), &end, 0.0, FALSE, 0.0, 0.0);
}

static void set_status(App *app, const char *msg) {
  gtk_label_set_text(GTK_LABEL(app->status_label), msg);
}

static gchar *get_text_from_view(GtkWidget *text_view) {
  GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  GtkTextIter start, end;
  gtk_text_buffer_get_start_iter(buf, &start);
  gtk_text_buffer_get_end_iter(buf, &end);
  return gtk_text_buffer_get_text(buf, &start, &end, FALSE);
}

static void set_text_to_view(GtkWidget *text_view, const gchar *text) {
  GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_view));
  gtk_text_buffer_set_text(buf, text ? text : "", -1);
}

static int selected_mode_value(App *app) {
  int idx = gtk_combo_box_get_active(GTK_COMBO_BOX(app->combo_mode));
  if (idx < 0 || idx >= (int)(sizeof(k_mode_values) / sizeof(k_mode_values[0]))) return 0;
  return k_mode_values[idx];
}

static void set_mode_combo(App *app, int mode) {
  for (int i = 0; i < (int)(sizeof(k_mode_values) / sizeof(k_mode_values[0])); i++) {
    if (k_mode_values[i] == mode) {
      gtk_combo_box_set_active(GTK_COMBO_BOX(app->combo_mode), i);
      return;
    }
  }
}

static gchar *upsert_param_line(const gchar *content, const gchar *key, const gchar *value) {
  GString *out = g_string_new("");
  gboolean found = FALSE;
  gchar **lines = g_strsplit(content ? content : "", "\n", -1);

  for (int i = 0; lines[i] != NULL; i++) {
    gchar *trim = g_strdup(lines[i]);
    g_strstrip(trim);
    if (g_str_has_prefix(trim, key)) {
      const gchar *p = trim + strlen(key);
      if (*p == '=' || g_str_has_prefix(p, " =")) {
        g_string_append_printf(out, "%s=%s\n", key, value);
        found = TRUE;
      } else {
        g_string_append(out, lines[i]);
        g_string_append_c(out, '\n');
      }
    } else {
      g_string_append(out, lines[i]);
      g_string_append_c(out, '\n');
    }
    g_free(trim);
  }

  if (!found) g_string_append_printf(out, "%s=%s\n", key, value);
  g_strfreev(lines);
  return g_string_free(out, FALSE);
}

static gchar *extract_param_value(const gchar *content, const gchar *key) {
  gchar **lines = g_strsplit(content ? content : "", "\n", -1);
  gchar *result = NULL;

  for (int i = 0; lines[i] != NULL; i++) {
    gchar *line = g_strdup(lines[i]);
    gchar *trim = g_strstrip(line);
    if (!g_str_has_prefix(trim, key)) {
      g_free(line);
      continue;
    }

    gchar *p = trim + strlen(key);
    while (*p == ' ' || *p == '\t') p++;
    if (*p != '=') {
      g_free(line);
      continue;
    }
    p++;
    while (*p == ' ' || *p == '\t') p++;
    gchar *comment = strstr(p, "//");
    if (comment) *comment = '\0';
    g_strstrip(p);
    result = g_strdup(p);
    g_free(line);
    break;
  }

  g_strfreev(lines);
  return result;
}

static gchar *build_inputpar_path(App *app) {
  const gchar *run_dir = gtk_entry_get_text(GTK_ENTRY(app->entry_run_dir));
  return g_build_filename(run_dir, "inputpar.txt", NULL);
}

static void sync_form_from_text(App *app) {
  gchar *content = get_text_from_view(app->text_input);
  for (int i = 0; i < app->nfields; i++) {
    gchar *v = extract_param_value(content, app->fields[i].key);
    if (v) {
      gtk_entry_set_text(GTK_ENTRY(app->fields[i].entry), v);
    } else if (g_strcmp0(app->fields[i].key, "eachopt") == 0) {
      gtk_entry_set_text(GTK_ENTRY(app->fields[i].entry), "0");
    } else {
      gtk_entry_set_text(GTK_ENTRY(app->fields[i].entry), "");
    }
    g_free(v);
  }
  gchar *mode = extract_param_value(content, "mode");
  if (mode && *mode) set_mode_combo(app, atoi(mode));
  g_free(mode);
  g_free(content);
}

static void sync_text_from_form(App *app) {
  gchar *content = get_text_from_view(app->text_input);
  gchar *next = content;

  for (int i = 0; i < app->nfields; i++) {
    const gchar *val = gtk_entry_get_text(GTK_ENTRY(app->fields[i].entry));
    if (!val || !*val) continue;
    gchar *patched = upsert_param_line(next, app->fields[i].key, val);
    g_free(next);
    next = patched;
  }

  gchar mode_str[32];
  g_snprintf(mode_str, sizeof(mode_str), "%d", selected_mode_value(app));
  gchar *patched_mode = upsert_param_line(next, "mode", mode_str);
  g_free(next);
  set_text_to_view(app->text_input, patched_mode);
  g_free(patched_mode);
}

static gboolean save_inputpar(App *app, GError **err) {
  sync_text_from_form(app);
  gchar *path = build_inputpar_path(app);
  gchar *content = get_text_from_view(app->text_input);
  gboolean ok = g_file_set_contents(path, content, -1, err);
  g_free(path);
  g_free(content);
  return ok;
}

static gboolean run_plot_hook(App *app) {
  const gchar *cmd = gtk_entry_get_text(GTK_ENTRY(app->entry_plot_cmd));
  const gchar *run_dir = gtk_entry_get_text(GTK_ENTRY(app->entry_run_dir));
  if (!cmd || !*cmd) return TRUE;

  gchar *argv[] = {(gchar *)"/bin/bash", (gchar *)"-lc", (gchar *)cmd, NULL};
  gchar *out = NULL;
  gchar *err = NULL;
  gint exit_status = 0;
  GError *gerr = NULL;

  append_log(app, "[GUI] Running plot hook...\n");
  if (!g_spawn_sync(run_dir, argv, NULL, 0, NULL, NULL, &out, &err, &exit_status, &gerr)) {
    gchar *msg = g_strdup_printf("[GUI] Plot hook failed to start: %s\n", gerr ? gerr->message : "unknown");
    append_log(app, msg);
    g_free(msg);
    if (gerr) g_error_free(gerr);
    g_free(out);
    g_free(err);
    return FALSE;
  }

  if (out && *out) append_log(app, out);
  if (err && *err) append_log(app, err);
  if (!WIFEXITED(exit_status) || WEXITSTATUS(exit_status) != 0) {
    append_log(app, "[GUI] Plot hook returned non-zero status.\n");
    g_free(out);
    g_free(err);
    return FALSE;
  }
  append_log(app, "[GUI] Plot hook completed.\n");
  g_free(out);
  g_free(err);
  return TRUE;
}

static void on_choose_path(GtkButton *btn, gpointer user_data) {
  PathPickerData *pp = (PathPickerData *)user_data;
  GtkFileChooserAction action = pp->select_folder ? GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER : GTK_FILE_CHOOSER_ACTION_OPEN;
  const gchar *title = pp->select_folder ? "Select Directory" : "Select File";

  GtkWidget *dialog = gtk_file_chooser_dialog_new(
      title, GTK_WINDOW(pp->window), action, "_Cancel", GTK_RESPONSE_CANCEL, "_Select", GTK_RESPONSE_ACCEPT, NULL);
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
    char *path = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
    gtk_entry_set_text(GTK_ENTRY(pp->entry), path);
    g_free(path);
  }
  gtk_widget_destroy(dialog);
  (void)btn;
}

static gboolean on_io_ready(GIOChannel *source, GIOCondition condition, gpointer data) {
  App *app = (App *)data;
  if (condition & (G_IO_HUP | G_IO_ERR | G_IO_NVAL)) return FALSE;

  gchar buf[2048];
  gsize bytes_read = 0;
  GError *err = NULL;
  GIOStatus st = g_io_channel_read_chars(source, buf, sizeof(buf) - 1, &bytes_read, &err);
  if (st == G_IO_STATUS_ERROR) {
    append_log(app, "[GUI] Failed to read process output.\n");
    if (err) g_error_free(err);
    return FALSE;
  }
  if (bytes_read > 0) {
    buf[bytes_read] = '\0';
    append_log(app, buf);
  }
  return TRUE;
}

static void cleanup_process_io(App *app) {
  if (app->stdout_watch_id) g_source_remove(app->stdout_watch_id);
  if (app->stderr_watch_id) g_source_remove(app->stderr_watch_id);
  app->stdout_watch_id = 0;
  app->stderr_watch_id = 0;
  if (app->stdout_ch) g_io_channel_unref(app->stdout_ch);
  if (app->stderr_ch) g_io_channel_unref(app->stderr_ch);
  app->stdout_ch = NULL;
  app->stderr_ch = NULL;
}

static void on_child_exit(GPid pid, gint status, gpointer data) {
  App *app = (App *)data;
  cleanup_process_io(app);
  g_spawn_close_pid(pid);
  app->child_pid = 0;
  app->running = FALSE;
  gtk_widget_set_sensitive(app->btn_run, TRUE);
  gtk_widget_set_sensitive(app->btn_stop, FALSE);

  if (WIFEXITED(status)) {
    gchar *msg = g_strdup_printf("[GUI] Process exited with code %d.\n", WEXITSTATUS(status));
    append_log(app, msg);
    g_free(msg);
  } else if (WIFSIGNALED(status)) {
    gchar *msg = g_strdup_printf("[GUI] Process terminated by signal %d.\n", WTERMSIG(status));
    append_log(app, msg);
    g_free(msg);
  } else {
    append_log(app, "[GUI] Process finished.\n");
  }
  set_status(app, "Idle");
}

static void on_load_input(GtkButton *btn, gpointer data) {
  App *app = (App *)data;
  gchar *path = build_inputpar_path(app);
  gchar *content = NULL;
  gsize len = 0;
  GError *err = NULL;

  if (!g_file_get_contents(path, &content, &len, &err)) {
    gchar *msg = g_strdup_printf("Failed to read %s: %s", path, err ? err->message : "unknown error");
    set_status(app, msg);
    append_log(app, "[GUI] ");
    append_log(app, msg);
    append_log(app, "\n");
    g_free(msg);
    if (err) g_error_free(err);
    g_free(path);
    return;
  }

  set_text_to_view(app->text_input, content);
  sync_form_from_text(app);
  set_status(app, "Loaded inputpar.txt");
  append_log(app, "[GUI] Loaded input parameters.\n");
  g_free(content);
  g_free(path);
  (void)btn;
}

static void on_save_input(GtkButton *btn, gpointer data) {
  App *app = (App *)data;
  GError *err = NULL;
  if (save_inputpar(app, &err)) {
    set_status(app, "Saved inputpar.txt");
    append_log(app, "[GUI] Saved input parameters.\n");
  } else {
    gchar *msg = g_strdup_printf("Failed to save inputpar.txt: %s", err ? err->message : "unknown error");
    set_status(app, msg);
    append_log(app, "[GUI] ");
    append_log(app, msg);
    append_log(app, "\n");
    g_free(msg);
  }
  if (err) g_error_free(err);
  (void)btn;
}

static void on_text_to_form(GtkButton *btn, gpointer data) {
  App *app = (App *)data;
  sync_form_from_text(app);
  set_status(app, "Updated form from text");
  (void)btn;
}

static void on_form_to_text(GtkButton *btn, gpointer data) {
  App *app = (App *)data;
  sync_text_from_form(app);
  set_status(app, "Applied form values to text");
  (void)btn;
}

static void on_plot_clicked(GtkButton *btn, gpointer data) {
  App *app = (App *)data;
  if (run_plot_hook(app)) set_status(app, "Plot command finished");
  else set_status(app, "Plot command failed");
  (void)btn;
}

static void on_stop_run(GtkButton *btn, gpointer data) {
  App *app = (App *)data;
  if (app->running && app->child_pid > 0) {
    kill(app->child_pid, SIGTERM);
    append_log(app, "[GUI] Sent SIGTERM to running process.\n");
    set_status(app, "Stopping...");
  }
  (void)btn;
}

static void on_start_run(GtkButton *btn, gpointer data) {
  App *app = (App *)data;
  if (app->running) return;

  const gchar *run_dir = gtk_entry_get_text(GTK_ENTRY(app->entry_run_dir));
  const gchar *bin = gtk_entry_get_text(GTK_ENTRY(app->entry_bin));
  int nproc = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(app->spin_nproc));
  GError *err = NULL;

  if (!g_file_test(run_dir, G_FILE_TEST_IS_DIR)) {
    set_status(app, "Run directory is invalid.");
    append_log(app, "[GUI] Invalid run directory.\n");
    return;
  }
  if (!g_file_test(bin, G_FILE_TEST_IS_REGULAR)) {
    set_status(app, "SMIwiz binary path is invalid.");
    append_log(app, "[GUI] Invalid SMIwiz binary.\n");
    return;
  }
  if (!save_inputpar(app, &err)) {
    gchar *msg = g_strdup_printf("Save failed before run: %s", err ? err->message : "unknown error");
    set_status(app, msg);
    append_log(app, "[GUI] ");
    append_log(app, msg);
    append_log(app, "\n");
    g_free(msg);
    if (err) g_error_free(err);
    return;
  }

  gchar *input_path = build_inputpar_path(app);
  gchar *cmd = g_strdup_printf("mpirun -n %d \"%s\" $(cat \"%s\")", nproc, bin, input_path);
  gchar *argv[] = {(gchar *)"/bin/bash", (gchar *)"-lc", cmd, NULL};
  gint out_fd = -1, err_fd = -1;

  append_log(app, "[GUI] Running: ");
  append_log(app, cmd);
  append_log(app, "\n");

  if (!g_spawn_async_with_pipes(run_dir, argv, NULL, G_SPAWN_DO_NOT_REAP_CHILD, NULL, NULL, &app->child_pid, NULL, &out_fd, &err_fd, &err)) {
    gchar *msg = g_strdup_printf("Failed to start process: %s", err ? err->message : "unknown error");
    set_status(app, msg);
    append_log(app, "[GUI] ");
    append_log(app, msg);
    append_log(app, "\n");
    g_free(msg);
    if (err) g_error_free(err);
    g_free(input_path);
    g_free(cmd);
    return;
  }

  app->stdout_ch = g_io_channel_unix_new(out_fd);
  app->stderr_ch = g_io_channel_unix_new(err_fd);
  g_io_channel_set_encoding(app->stdout_ch, NULL, NULL);
  g_io_channel_set_encoding(app->stderr_ch, NULL, NULL);
  g_io_channel_set_flags(app->stdout_ch, G_IO_FLAG_NONBLOCK, NULL);
  g_io_channel_set_flags(app->stderr_ch, G_IO_FLAG_NONBLOCK, NULL);
  app->stdout_watch_id = g_io_add_watch(app->stdout_ch, G_IO_IN | G_IO_HUP | G_IO_ERR, on_io_ready, app);
  app->stderr_watch_id = g_io_add_watch(app->stderr_ch, G_IO_IN | G_IO_HUP | G_IO_ERR, on_io_ready, app);
  g_child_watch_add(app->child_pid, on_child_exit, app);

  app->running = TRUE;
  gtk_widget_set_sensitive(app->btn_run, FALSE);
  gtk_widget_set_sensitive(app->btn_stop, TRUE);
  set_status(app, "Running");

  g_free(input_path);
  g_free(cmd);
  (void)btn;
}

static GtkWidget *make_labeled_entry(GtkWidget *grid, const char *label, int row, GtkWidget **entry_out) {
  GtkWidget *lbl = gtk_label_new(label);
  gtk_widget_set_halign(lbl, GTK_ALIGN_START);
  GtkWidget *entry = gtk_entry_new();
  gtk_grid_attach(GTK_GRID(grid), lbl, 0, row, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), entry, 1, row, 1, 1);
  *entry_out = entry;
  return entry;
}

static void add_form_field(App *app, GtkWidget *grid, int idx, int slot, const char *key, const char *label) {
  int col_group = slot % 3;
  int row = slot / 3;
  int col = col_group * 2;
  GtkWidget *lbl = gtk_label_new(label);
  GtkWidget *entry = gtk_entry_new();
  gtk_widget_set_halign(lbl, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), lbl, col, row, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), entry, col + 1, row, 1, 1);
  app->fields[idx].key = key;
  app->fields[idx].label = label;
  app->fields[idx].entry = entry;
}

static void populate_default_paths(App *app) {
  gchar *cwd = g_get_current_dir();
  gchar *project = NULL;
  gchar *base = g_path_get_basename(cwd);
  if (g_strcmp0(base, "GUI") == 0) project = g_path_get_dirname(cwd);
  else project = g_strdup(cwd);

  gchar *bin = g_build_filename(project, "bin", "SMIwiz", NULL);
  gchar *run_dir = g_build_filename(project, "run_fwd", NULL);
  gtk_entry_set_text(GTK_ENTRY(app->entry_project), project);
  gtk_entry_set_text(GTK_ENTRY(app->entry_bin), bin);
  gtk_entry_set_text(GTK_ENTRY(app->entry_run_dir), run_dir);

  g_free(run_dir);
  g_free(bin);
  g_free(base);
  g_free(project);
  g_free(cwd);
}

int main(int argc, char **argv) {
  gtk_init(&argc, &argv);

  App app = {0};
  app.nfields = 16;

  app.window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(app.window), "SMIwiz GTK Launcher");
  gtk_window_set_default_size(GTK_WINDOW(app.window), 1300, 860);
  g_signal_connect(app.window, "destroy", G_CALLBACK(gtk_main_quit), NULL);

  GtkWidget *root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_container_set_border_width(GTK_CONTAINER(root), 8);
  gtk_container_add(GTK_CONTAINER(app.window), root);

  GtkWidget *main_paned = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
  gtk_box_pack_start(GTK_BOX(root), main_paned, TRUE, TRUE, 0);

  GtkWidget *left = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_paned_add1(GTK_PANED(main_paned), left);

  GtkWidget *right = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
  gtk_paned_add2(GTK_PANED(main_paned), right);

  GtkWidget *grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(grid), 6);
  gtk_grid_set_column_spacing(GTK_GRID(grid), 8);
  gtk_box_pack_start(GTK_BOX(left), grid, FALSE, FALSE, 0);
  make_labeled_entry(grid, "Project Dir", 0, &app.entry_project);
  make_labeled_entry(grid, "Run Dir", 1, &app.entry_run_dir);
  make_labeled_entry(grid, "SMIwiz Binary", 2, &app.entry_bin);

  GtkWidget *btn_project = gtk_button_new_with_label("Browse");
  GtkWidget *btn_run_dir = gtk_button_new_with_label("Browse");
  GtkWidget *btn_bin = gtk_button_new_with_label("Browse");
  gtk_grid_attach(GTK_GRID(grid), btn_project, 2, 0, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), btn_run_dir, 2, 1, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), btn_bin, 2, 2, 1, 1);

  GtkWidget *lbl_np = gtk_label_new("MPI Ranks");
  app.spin_nproc = gtk_spin_button_new_with_range(1, 1024, 1);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(app.spin_nproc), 2);
  gtk_widget_set_halign(lbl_np, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), lbl_np, 0, 3, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), app.spin_nproc, 1, 3, 1, 1);

  GtkWidget *lbl_mode = gtk_label_new("Mode");
  app.combo_mode = gtk_combo_box_text_new();
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "0 Forward modeling");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "1 FWI");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "2 RTM");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "3 Data-domain LSRTM");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "4 FWI gradient");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "5 Source inversion");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "6 ADCIG");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "7 PSF Hessian");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "8 Migration deconvolution (PCG)");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "9 Migration deconvolution (FFT-Wiener)");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_mode), "10 Up-down wavefield separation");
  gtk_combo_box_set_active(GTK_COMBO_BOX(app.combo_mode), 0);
  gtk_widget_set_halign(lbl_mode, GTK_ALIGN_START);
  gtk_grid_attach(GTK_GRID(grid), lbl_mode, 0, 4, 1, 1);
  gtk_grid_attach(GTK_GRID(grid), app.combo_mode, 1, 4, 2, 1);

  GtkWidget *form_frame = gtk_frame_new("Common Parameters (Form)");
  gtk_box_pack_start(GTK_BOX(left), form_frame, FALSE, FALSE, 0);
  GtkWidget *form_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
  gtk_container_set_border_width(GTK_CONTAINER(form_box), 6);
  gtk_container_add(GTK_CONTAINER(form_frame), form_box);

  GtkWidget *file_frame = gtk_frame_new("File Names");
  gtk_box_pack_start(GTK_BOX(form_box), file_frame, FALSE, FALSE, 0);
  GtkWidget *file_grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(file_grid), 4);
  gtk_grid_set_column_spacing(GTK_GRID(file_grid), 8);
  gtk_container_set_border_width(GTK_CONTAINER(file_grid), 6);
  gtk_container_add(GTK_CONTAINER(file_frame), file_grid);
  add_form_field(&app, file_grid, 0, 0, "acquifile", "acquifile");
  add_form_field(&app, file_grid, 1, 1, "vpfile", "vpfile");
  add_form_field(&app, file_grid, 2, 2, "rhofile", "rhofile");
  add_form_field(&app, file_grid, 3, 3, "stffile", "stffile");

  GtkWidget *other_frame = gtk_frame_new("Other Parameters");
  gtk_box_pack_start(GTK_BOX(form_box), other_frame, FALSE, FALSE, 0);
  GtkWidget *other_grid = gtk_grid_new();
  gtk_grid_set_row_spacing(GTK_GRID(other_grid), 4);
  gtk_grid_set_column_spacing(GTK_GRID(other_grid), 8);
  gtk_container_set_border_width(GTK_CONTAINER(other_grid), 6);
  gtk_container_add(GTK_CONTAINER(other_frame), other_grid);
  add_form_field(&app, other_grid, 15, 0, "eachopt", "eachopt");
  add_form_field(&app, other_grid, 4, 1, "freesurf", "freesurf");
  add_form_field(&app, other_grid, 5, 2, "nt", "nt");
  add_form_field(&app, other_grid, 6, 3, "dt", "dt");
  add_form_field(&app, other_grid, 7, 4, "nb", "nb");
  add_form_field(&app, other_grid, 8, 5, "n1", "n1");
  add_form_field(&app, other_grid, 9, 6, "n2", "n2");
  add_form_field(&app, other_grid, 10, 7, "n3", "n3");
  add_form_field(&app, other_grid, 11, 8, "d1", "d1");
  add_form_field(&app, other_grid, 12, 9, "d2", "d2");
  add_form_field(&app, other_grid, 13, 10, "d3", "d3");
  add_form_field(&app, other_grid, 14, 11, "order", "order");

  GtkWidget *sync_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
  GtkWidget *btn_load = gtk_button_new_with_label("Load inputpar.txt");
  GtkWidget *btn_save = gtk_button_new_with_label("Save inputpar.txt");
  GtkWidget *btn_text_to_form = gtk_button_new_with_label("Text -> Form");
  GtkWidget *btn_form_to_text = gtk_button_new_with_label("Form -> Text");
  gtk_box_pack_start(GTK_BOX(sync_row), btn_load, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(sync_row), btn_save, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(sync_row), btn_text_to_form, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(sync_row), btn_form_to_text, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(form_box), sync_row, FALSE, FALSE, 0);

  GtkWidget *run_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
  app.btn_run = gtk_button_new_with_label("Run");
  app.btn_stop = gtk_button_new_with_label("Stop");
  gtk_widget_set_sensitive(app.btn_stop, FALSE);
  gtk_box_pack_start(GTK_BOX(run_row), app.btn_run, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(run_row), app.btn_stop, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(right), run_row, FALSE, FALSE, 0);

  GtkWidget *plot_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
  GtkWidget *btn_plot = gtk_button_new_with_label("Plot");
  app.entry_plot_cmd = gtk_entry_new();
  gtk_entry_set_text(GTK_ENTRY(app.entry_plot_cmd), "bash plot.sh");
  GtkWidget *frame_input = gtk_frame_new("inputpar.txt");
  GtkWidget *sw_input = gtk_scrolled_window_new(NULL, NULL);
  app.text_input = gtk_text_view_new();
  gtk_container_add(GTK_CONTAINER(sw_input), app.text_input);
  gtk_container_add(GTK_CONTAINER(frame_input), sw_input);
  gtk_box_pack_start(GTK_BOX(left), frame_input, TRUE, TRUE, 0);
  GtkWidget *frame_log = gtk_frame_new("Run Log");
  GtkWidget *sw_log = gtk_scrolled_window_new(NULL, NULL);
  app.text_log = gtk_text_view_new();
  gtk_text_view_set_editable(GTK_TEXT_VIEW(app.text_log), FALSE);
  gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(app.text_log), FALSE);
  gtk_container_add(GTK_CONTAINER(sw_log), app.text_log);
  gtk_container_add(GTK_CONTAINER(frame_log), sw_log);
  gtk_box_pack_start(GTK_BOX(right), frame_log, TRUE, TRUE, 0);

  gtk_box_pack_start(GTK_BOX(plot_row), btn_plot, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(plot_row), gtk_label_new("Plot Command"), FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(plot_row), app.entry_plot_cmd, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(right), plot_row, FALSE, FALSE, 0);

  gtk_paned_set_position(GTK_PANED(main_paned), 880);

  app.status_label = gtk_label_new("Idle");
  gtk_widget_set_halign(app.status_label, GTK_ALIGN_START);
  gtk_box_pack_start(GTK_BOX(root), app.status_label, FALSE, FALSE, 0);

  populate_default_paths(&app);
  on_load_input(GTK_BUTTON(btn_load), &app);

  PathPickerData *project_choose = g_new0(PathPickerData, 1);
  PathPickerData *run_choose = g_new0(PathPickerData, 1);
  PathPickerData *bin_choose = g_new0(PathPickerData, 1);
  project_choose->entry = app.entry_project;
  project_choose->window = app.window;
  project_choose->select_folder = TRUE;
  run_choose->entry = app.entry_run_dir;
  run_choose->window = app.window;
  run_choose->select_folder = TRUE;
  bin_choose->entry = app.entry_bin;
  bin_choose->window = app.window;
  bin_choose->select_folder = FALSE;

  g_signal_connect(btn_project, "clicked", G_CALLBACK(on_choose_path), project_choose);
  g_signal_connect(btn_run_dir, "clicked", G_CALLBACK(on_choose_path), run_choose);
  g_signal_connect(btn_bin, "clicked", G_CALLBACK(on_choose_path), bin_choose);
  g_signal_connect(btn_load, "clicked", G_CALLBACK(on_load_input), &app);
  g_signal_connect(btn_save, "clicked", G_CALLBACK(on_save_input), &app);
  g_signal_connect(btn_text_to_form, "clicked", G_CALLBACK(on_text_to_form), &app);
  g_signal_connect(btn_form_to_text, "clicked", G_CALLBACK(on_form_to_text), &app);
  g_signal_connect(btn_plot, "clicked", G_CALLBACK(on_plot_clicked), &app);
  g_signal_connect(app.btn_run, "clicked", G_CALLBACK(on_start_run), &app);
  g_signal_connect(app.btn_stop, "clicked", G_CALLBACK(on_stop_run), &app);

  gtk_widget_show_all(app.window);
  gtk_main();
  cleanup_process_io(&app);
  return 0;
}
